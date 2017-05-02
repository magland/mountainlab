/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/4/2016
*******************************************************/

#include "mdaconvert.h"
#include "mdaio.h"

#include <QFile>
#include <QFileInfo>
#include <QTime>
#include <mda.h>
#include <stdio.h>

#include "mlcommon.h"

bigint get_mda_dtype(QString format);

struct working_data {
    ~working_data()
    {
        if (inf)
            fclose(inf);
        if (outf)
            fclose(outf);
        if (buf)
            free(buf);
    }

    MDAIO_HEADER HH_in;
    MDAIO_HEADER HH_out;
    FILE* inf = 0;
    FILE* outf = 0;
    void* buf = 0;
};

bool copy_data(working_data& D, bigint N);
bigint get_num_bytes_per_entry(bigint dtype);
bool convert_ncs(const mdaconvert_opts& opts);
bool convert_nrd(const mdaconvert_opts& opts);

bool mdaconvert(const mdaconvert_opts& opts_in)
{
    mdaconvert_opts opts = opts_in;

    //default inputs in case input format is mda, csv, or dat
    if (opts.input_format == "mda") {
        if (!opts.input_dtype.isEmpty()) {
            qWarning() << "input-dtype should not be specified for input type mda";
            return false;
        }
        if (!opts.dims.isEmpty()) {
            qWarning() << "dims should not be specified for input type mda";
            return false;
        }
    }
    else if (opts.input_format == "csv") {
        if (!opts.dims.isEmpty()) {
            qWarning() << "dims should not be specified for input type csv";
            return false;
        }
    }
    else if (opts.input_format == "nrd") {
        if (!opts.dims.isEmpty()) {
            qWarning() << "dims should not be specified for input type nrd";
            return false;
        }
        if (opts.num_channels == 0) {
            qWarning() << "Number of channels must be specified for type nrd. Use, for example, --num_channels=32";
            return false;
        }
    }
    else if (opts.input_format == "ncs") {
        if (!opts.dims.isEmpty()) {
            qWarning() << "dims should not be specified for input type ncs";
            return false;
        }
    }
    else if (opts.input_format == "raw_timeseries") {
        if (opts.input_dtype.isEmpty()) {
            qWarning() << "Input datatype (e.g., --dtype=int16)  must be specified";
            return false;
        }
        if (!opts.dims.isEmpty()) {
            qWarning() << "dims should not be specified for input type raw_timeseries. Use, e.g., --num_channels=32 instead";
            return false;
        }
        if (opts.num_channels == 0) {
            qWarning() << "Number of channels must be specified for type raw_timeeries. Use, for example, --num_channels=32";
            return false;
        }
        opts.dims.clear();
        opts.dims << opts.num_channels;
        bigint input_size = QFileInfo(opts.input_path).size();
        bigint dtype_size = get_num_bytes_per_entry(get_mda_dtype(opts.input_dtype));
        bigint num_timepoints = input_size / (dtype_size * opts.num_channels);
        if (input_size != num_timepoints * dtype_size * opts.num_channels) {
            qWarning() << QString("The size of the file (%1) is not compatible with this datatype (%2) and number of channels (%3).").arg(input_size).arg(opts.input_dtype).arg(opts.num_channels);
            return false;
        }
        printf("Auto-setting number of timepoints to %ld\n", num_timepoints);
        opts.dims << num_timepoints;
    }
    else {
        if (opts.input_dtype.isEmpty()) {
            qWarning() << "Input datatype (e.g., --dtype=int16) must be specified";
            return false;
        }
        if (opts.dims.isEmpty()) {
            qWarning() << "Dimensions must be specified";
            return false;
        }
    }

    //check file existence, etc
    if (!QFile::exists(opts.input_path)) {
        qWarning() << "Input file does not exist: " + opts.input_path;
        return false;
    }
    if (QFile::exists(opts.output_path)) {
        if (!QFile::remove(opts.output_path)) {
            qWarning() << "Unable to remove output file: " + opts.output_path;
            return false;
        }
    }

    if (opts.input_format == "csv") {
        if (opts.output_format != "mda") {
            qWarning() << "Invalid output type for input type csv";
            return false;
        }
        QString txt = TextFile::read(opts.input_path);
        QStringList lines = txt.split("\n");
        lines = lines.mid(opts.input_num_header_rows);
        QString line0 = lines.value(0);
        bigint num_cols = line0.split(",").mid(opts.input_num_header_cols).count();
        if (lines.value(lines.count() - 1).trimmed().isEmpty())
            lines = lines.mid(0, lines.count() - 1);
        bigint num_rows = lines.count();
        Mda out(num_cols, num_rows);
        for (bigint r = 0; r < lines.count(); r++) {
            QString line = lines[r];
            QStringList vals = line.split(",").mid(opts.input_num_header_cols);
            QList<double> vals2;
            for (bigint j = 0; j < vals.count(); j++) {
                vals2 << vals[j].trimmed().toDouble();
            }
            for (bigint c = 0; c < num_cols; c++) {
                out.setValue(vals2.value(c), c, r);
            }
        }
        out.write32(opts.output_path);
        return true;
    }
    if (opts.input_format == "nrd") {
        return convert_nrd(opts);
    }
    if (opts.input_format == "ncs") {
        return convert_ncs(opts);
    }

    // initialize working data
    working_data D;
    for (bigint i = 0; i < MDAIO_MAX_DIMS; i++) {
        D.HH_in.dims[i] = 1;
        D.HH_out.dims[i] = 1;
    }

    // open input file, read header if mda
    D.inf = fopen(opts.input_path.toLatin1().data(), "rb");
    if (!D.inf) {
        qWarning() << "Unable to open input file for reading: " + opts.input_path;
        return false;
    }
    // read header for format = mda
    if (opts.input_format == "mda") {
        if (!mda_read_header(&D.HH_in, D.inf)) {
            qWarning() << "Error reading input header";
            return false;
        }
    }
    else {
        D.HH_in.data_type = get_mda_dtype(opts.input_dtype);
        D.HH_in.header_size = 0;
        D.HH_in.num_bytes_per_entry = get_num_bytes_per_entry(D.HH_in.data_type);
        if (!D.HH_in.num_bytes_per_entry) {
            qWarning() << "Unable to determine input number of bytes per entry.";
            return false;
        }
        if (opts.dims.value(opts.dims.count() - 1) == 0) {
            qWarning() << "Final dimension of 0 is no longer supported because there was apparently a bug. If your array is two-dimensional you may want to use --input_format=raw_timeseries.\n";
            return false;
            /*
            //auto-determine final dimension
            bigint size0 = QFileInfo(opts.input_path).size();
            bigint tmp = size0 / D.HH_in.num_bytes_per_entry;
            for (bigint j = 0; j < opts.dims.count() - 1; j++) {
                tmp /= opts.dims[j];
            }
            printf("Auto-setting final dimension to %d\n", tmp);
            opts.dims[opts.dims.count() - 1] = tmp;
            */
        }
        D.HH_in.num_dims = opts.dims.count();
        for (bigint i = 0; i < D.HH_in.num_dims; i++) {
            D.HH_in.dims[i] = opts.dims[i];
        }
    }

    if (D.HH_in.num_dims < 2) {
        qWarning() << "Number of dimensions must be at least 2.";
        return false;
    }
    if (D.HH_in.num_dims > 6) {
        qWarning() << "Number of dimensions can be at most 6.";
        return false;
    }

    // output header
    D.HH_out.data_type = D.HH_in.data_type;
    if (!opts.output_dtype.isEmpty()) {
        D.HH_out.data_type = get_mda_dtype(opts.output_dtype);
        if (!D.HH_out.data_type) {
            qWarning() << "Invalid output datatype:" << opts.output_dtype;
            return false;
        }
    }
    D.HH_out.num_bytes_per_entry = get_num_bytes_per_entry(D.HH_out.data_type);
    D.HH_out.num_dims = D.HH_in.num_dims;
    bigint dim_prod = 1;
    for (bigint i = 0; i < D.HH_out.num_dims; i++) {
        D.HH_out.dims[i] = D.HH_in.dims[i];
        dim_prod *= D.HH_in.dims[i];
    }

    //open output file, write header if mda
    D.outf = fopen(opts.output_path.toLatin1().data(), "wb");
    if (!D.outf) {
        qWarning() << "Unable to open output file for writing: " + opts.output_path;
        return false;
    }
    if (opts.output_format == "mda") {
        if (!mda_write_header(&D.HH_out, D.outf)) {
            qWarning() << "Error writing output header";
            return false;
        }
    }

    //check input file size
    if (opts.check_input_file_size) {
        bigint expected_input_file_size = D.HH_in.num_bytes_per_entry * dim_prod + D.HH_in.header_size;
        bigint actual_input_file_size = QFileInfo(opts.input_path).size();
        if (actual_input_file_size != expected_input_file_size) {
            qWarning() << QString("Unexpected input file size: Expected/Actual=%1/%2 bytes").arg(expected_input_file_size).arg(actual_input_file_size);
            qWarning() << "Try using --allow-subset";
            return false;
        }
    }

    bool ret = true;
    bigint chunk_size = 1e7;
    QTime timer;
    timer.start();
    for (bigint ii = 0; ii < dim_prod; ii += chunk_size) {
        if (timer.elapsed() > 1000) {
            bigint pct = (bigint)((ii * 1.0 / dim_prod) * 100);
            printf("%ld%%\n", pct);
            timer.restart();
        }
        bigint NN = qMin(chunk_size, dim_prod - ii);
        if (!copy_data(D, NN)) {
            ret = false;
            break;
        }
    }

    if (!ret) {
        fclose(D.outf);
        D.outf = 0;
        QFile::remove(opts.output_path);
    }
    return false;
}

bool copy_data(working_data& D, bigint N)
{
    if (!N)
        return true;
    if (D.buf) {
        free(D.buf);
        D.buf = 0;
    }
    bigint buf_size = N * D.HH_out.num_bytes_per_entry;
    D.buf = malloc(buf_size);
    if (!D.buf) {
        qWarning() << QString("Error in malloc of size %1 ").arg(buf_size);
        return false;
    }

    void* buf = D.buf;
    bigint num_read = 0;
    if (D.HH_out.data_type == MDAIO_TYPE_BYTE) {
        num_read = mda_read_byte((unsigned char*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_UINT16) {
        num_read = mda_read_uint16((quint16*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_INT16) {
        num_read = mda_read_int16((qint16*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_INT32) {
        num_read = mda_read_int32((qint32*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_FLOAT32) {
        num_read = mda_read_float32((float*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_FLOAT64) {
        num_read = mda_read_float64((double*)buf, &D.HH_in, N, D.inf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_UINT32) {
        num_read = mda_read_uint32((quint32*)buf, &D.HH_in, N, D.inf);
    }
    else {
        qWarning() << "Unsupported format for reading";
        return false;
    }

    if (num_read != N) {
        qWarning() << "Error reading data" << num_read << N;
        return false;
    }

    bigint num_written = 0;
    if (D.HH_out.data_type == MDAIO_TYPE_BYTE) {
        num_written = mda_write_byte((unsigned char*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_UINT16) {
        num_written = mda_write_uint16((quint16*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_INT16) {
        num_written = mda_write_int16((qint16*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_INT32) {
        num_written = mda_write_int32((qint32*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_FLOAT32) {
        num_written = mda_write_float32((float*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_FLOAT64) {
        num_written = mda_write_float64((double*)buf, &D.HH_out, N, D.outf);
    }
    else if (D.HH_out.data_type == MDAIO_TYPE_UINT32) {
        num_written = mda_write_uint32((quint32*)buf, &D.HH_out, N, D.outf);
    }
    else {
        qWarning() << "Unsupported format for writing";
        return false;
    }

    if (num_written != N) {
        qWarning() << "Error writing data" << num_written << N;
        return false;
    }

    return true;
}

/*
function test_read_ncs

% Open the file
F=fopen('CSC1.ncs');

% Skip the header
fseek(F,16*1024,'cof');

% Read all of the data
all_data=fread(F,inf,'int16');

% Compute the number of records
record_size_bytes=20+512*2;
data_size_bytes=length(all_data)*2;
num_records=floor(data_size_bytes/record_size_bytes);

% Use a reshape trick to remove the record headers (20 bytes each)
all_data=reshape(all_data(1:record_size_bytes/2*num_records),record_size_bytes/2,num_records);
all_data=all_data(11:end,:);
all_data=all_data(:);

% Close the file
fclose(F);

% Plot the first bit
figure; plot(all_data(1:2000));
*/

bool convert_ncs(const mdaconvert_opts& opts)
{
    FILE* inf = fopen(opts.input_path.toUtf8().data(), "r");
    if (!inf) {
        qWarning() << "Unable to open input file for reading: " + opts.input_path;
        return false;
    }
    FILE* outf = fopen(opts.output_path.toUtf8().data(), "w");
    if (!outf) {
        qWarning() << "Unable to open output file for writing: " + opts.output_path;
        fclose(inf);
        return false;
    }
    bigint size_of_input = QFileInfo(opts.input_path).size();
    bigint header_size = 16 * 1024;
    bigint record_size = 20 + 512 * 2;
    bigint num_records = (size_of_input - header_size) / record_size;
    fseek(inf, header_size, SEEK_SET); //skip the header
    MDAIO_HEADER H_out;
    H_out.data_type = MDAIO_TYPE_INT16;
    H_out.num_bytes_per_entry = 2;
    H_out.dims[0] = 1;
    H_out.dims[1] = num_records * 512;
    H_out.num_dims = 2;
    mda_write_header(&H_out, outf);
    char buffer[record_size];
    qint16* record_data = (qint16*)(&buffer[record_size - 512 * 2]);
    for (bigint i = 0; i < num_records; i++) {
        {
            bigint ret = fread(buffer, sizeof(char), record_size, inf);
            if (ret != record_size) {
                qWarning() << QString("Problem reading record %1 of %2 from .ncs file: ").arg(i).arg(num_records) + opts.input_path;
                fclose(inf);
                fclose(outf);
                return false;
            }
        }
        {
            bigint ret = mda_write_int16(record_data, &H_out, 512, outf);
            if (ret != 512) {
                qWarning() << QString("Problem writing record %1 of %2 from .ncs file: ").arg(i).arg(num_records) + opts.input_path;
                fclose(inf);
                fclose(outf);
                return false;
            }
        }
    }
    bigint test_num_bytes = fread(buffer, sizeof(char), record_size, inf);
    if (test_num_bytes != 0) {
        qWarning() << "Warning: extra bytes found in input file: " + opts.input_path;
    }

    fclose(inf);
    fclose(outf);

    return true;
}

/*

function test_read_nrd

fname='/home/spike/lilia/2017-01-25_12-51-32/DigitalLynxSXRawDataFile.nrd';
num_channels=96;

% Open the file
F=fopen(fname);

% Skip the header
fseek(F,16*1024,'cof');

% Read all of the data
all_data=fread(F,inf,'int32');

% Compute the number of records
record_size_bytes=32+40+4*num_channels;
data_size_bytes=length(all_data)*4;
num_records=floor(data_size_bytes/record_size_bytes);

% Use a reshape trick to extract the data part
all_data=reshape(all_data(1:record_size_bytes/4*num_records),record_size_bytes/4,num_records);
all_data=all_data(28+40+1:end-4,:);

% Close the file
fclose(F);

% Plot the first bit on the first bunch of channels
figure; plot(all_data(1:25,1:50000)');

*/

bool convert_nrd(const mdaconvert_opts& opts)
{
    FILE* inf = fopen(opts.input_path.toUtf8().data(), "r");
    if (!inf) {
        qWarning() << "Unable to open input file for reading: " + opts.input_path;
        return false;
    }
    FILE* outf = fopen(opts.output_path.toUtf8().data(), "w");
    if (!outf) {
        qWarning() << "Unable to open output file for writing: " + opts.output_path;
        fclose(inf);
        return false;
    }
    bigint size_of_input = QFileInfo(opts.input_path).size();
    bigint header_size = 16 * 1024;
    bigint record_size = 32 + 40 + 4 * opts.num_channels; //TODO: this should be auto-detected
    bigint num_records = (size_of_input - header_size) / record_size;
    fseek(inf, header_size, SEEK_SET); //skip the header
    MDAIO_HEADER H_out;
    H_out.data_type = MDAIO_TYPE_INT32;
    H_out.num_bytes_per_entry = 2;
    H_out.dims[0] = opts.num_channels;
    H_out.dims[1] = num_records;
    H_out.num_dims = 2;
    mda_write_header(&H_out, outf);
    char buffer[record_size];
    qint32* record_data = (qint32*)(&buffer[record_size - opts.num_channels * 4 - 4]);
    for (bigint i = 0; i < num_records; i++) {
        {
            bigint ret = fread(buffer, sizeof(char), record_size, inf);
            if (ret != record_size) {
                qWarning() << QString("Problem reading record %1 of %2 from .ncs file: ").arg(i).arg(num_records) + opts.input_path;
                fclose(inf);
                fclose(outf);
                return false;
            }
        }
        {
            bigint ret = mda_write_int32(record_data, &H_out, 512, outf);
            if (ret != 512) {
                qWarning() << QString("Problem writing record %1 of %2 from .ncs file: ").arg(i).arg(num_records) + opts.input_path;
                fclose(inf);
                fclose(outf);
                return false;
            }
        }
    }
    bigint test_num_bytes = fread(buffer, sizeof(char), record_size, inf);
    if (test_num_bytes != 0) {
        qWarning() << "Warning: extra bytes found in input file: " + opts.input_path;
    }

    fclose(inf);
    fclose(outf);

    return true;
}

bigint get_mda_dtype(QString format)
{
    if (format == "byte")
        return MDAIO_TYPE_BYTE;
    if (format == "int8")
        return MDAIO_TYPE_BYTE;
    if (format == "int16")
        return MDAIO_TYPE_INT16;
    if (format == "int32")
        return MDAIO_TYPE_INT32;
    if (format == "uint16")
        return MDAIO_TYPE_UINT16;
    if (format == "float32")
        return MDAIO_TYPE_FLOAT32;
    if (format == "float64")
        return MDAIO_TYPE_FLOAT64;
    if (format == "int32")
        return MDAIO_TYPE_INT32;
    return 0;
}

bigint get_num_bytes_per_entry(bigint dtype)
{
    switch (dtype) {
    case MDAIO_TYPE_BYTE:
        return 1;
        break;
    case MDAIO_TYPE_UINT16:
        return 2;
        break;
    case MDAIO_TYPE_INT16:
        return 2;
        break;
    case MDAIO_TYPE_INT32:
        return 4;
        break;
    case MDAIO_TYPE_FLOAT32:
        return 4;
        break;
    case MDAIO_TYPE_FLOAT64:
        return 8;
        break;
    case MDAIO_TYPE_UINT32:
        return 4;
        break;
    }

    return 0;
}
