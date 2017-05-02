/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/4/2016
*******************************************************/

#include "mdaconvert.h"
#include "mlcommon.h"

#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <diskreadmda.h>

void print_usage();

QString get_default_format(QString path);
QString get_info_string(const DiskReadMda& X);
QString get_json_header_string(const DiskReadMda& X);
bool extract_time_chunk(QString input_fname, QString output_fname, const QMap<QString, QVariant>& params);
bool extract_channels(QString input_fname, QString output_fname, const QMap<QString, QVariant>& params);

/// TODO, auto-calculate the last dimension

int main(int argc, char* argv[])
{
    CLParams params(argc, argv);

    if (params.unnamed_parameters.count() < 1) {
        print_usage();
        return 0;
    }

    QString arg1 = params.unnamed_parameters.value(0);
    QString arg2 = params.unnamed_parameters.value(1);
    QString arg3 = params.unnamed_parameters.value(2);

    if (arg1 == "extract_time_chunk") {
        if (extract_time_chunk(arg2, arg3, params.named_parameters))
            return 0;
        else
            return -1;
    }
    else if (arg1 == "extract_channels") {
        if (extract_channels(arg2, arg3, params.named_parameters))
            return 0;
        else
            return -1;
    }

    mdaconvert_opts opts;

    opts.input_path = params.unnamed_parameters.value(0);
    opts.input_dtype = params.named_parameters.value("input_dtype").toString();
    opts.input_format = params.named_parameters.value("input_format").toString();

    opts.output_path = params.unnamed_parameters.value(1);
    opts.output_dtype = params.named_parameters.value("output_dtype").toString();
    opts.output_format = params.named_parameters.value("output_format").toString();

    if (params.named_parameters.contains("allow-subset"))
        opts.check_input_file_size = false;

    if (opts.output_path.isEmpty()) {
        //if (opts.input_path.endsWith(".mda")) {
        DiskReadMda X(opts.input_path);
        if (params.named_parameters.contains("readheader")) {
            QString str = get_json_header_string(X);
            printf("%s\n", str.toUtf8().data());
            return 0;
        }
        else {
            QString str = get_info_string(X);
            printf("%s\n", str.toUtf8().data());
            return 0;
        }
        //}
        //else {
        //    print_usage();
        //    return -1;
        //}
    }

    if (params.named_parameters.contains("dtype")) {
        opts.input_dtype = params.named_parameters.value("dtype").toString();
        opts.output_dtype = params.named_parameters.value("dtype").toString();
    }

    if (params.named_parameters.contains("input_num_header_rows")) {
        opts.input_num_header_rows = params.named_parameters["input_num_header_rows"].toInt();
    }
    if (params.named_parameters.contains("input_num_header_cols")) {
        opts.input_num_header_cols = params.named_parameters["input_num_header_cols"].toInt();
    }

    if (opts.input_format.isEmpty()) {
        opts.input_format = get_default_format(opts.input_path);
    }
    if (opts.output_format.isEmpty()) {
        opts.output_format = get_default_format(opts.output_path);
    }

    QStringList dims_strlist = params.named_parameters.value("dims", "").toString().split("x", QString::SkipEmptyParts);
    foreach (QString str, dims_strlist) {
        opts.dims << str.toLong();
    }
    opts.num_channels = params.named_parameters.value("num_channels", 0).toInt();

    printf("Converting %s --> %s\n", opts.input_path.toLatin1().data(), opts.output_path.toLatin1().data());
    if (!mdaconvert(opts)) {
        return -1;
    }

    return 0;
}

QString get_suffix(QString path)
{
    int index1 = path.lastIndexOf(".");
    if (index1 < 0)
        return "";
    int index2 = path.lastIndexOf("/");
    if ((index2 >= 0) && (index2 > index1))
        return "";
    return path.mid(index1 + 1);
}

QString get_default_format(QString path)
{
    QString suf = get_suffix(path);
    if (suf == "mda")
        return "mda";
    else if (suf == "csv")
        return "csv";
    else if (suf.toLower() == "nrd") {
        return "nrd";
    }
    else if (suf.toLower() == "ncs") {
        return "ncs";
    }
    else
        return "raw";
}

void print_usage()
{
    printf("Example usages for converting between raw and mda formats:\n");
    printf("mdaconvert input.mda\n");
    printf("mdaconvert input.dat output.mda --dtype=uint16 --dims=32x100x44\n");
    printf("mdaconvert input.dat output.mda --dtype=int16 --input_format=raw_timeseries --num_channels=32\n");
    printf("mdaconvert input.file output.file --input_format=dat --input_dtype=float64 --output_format=mda --output_dtype=float32\n");
    printf("mdaconvert input.csv output.mda --input_num_header_rows=1 --input_num_header_cols=0\n");
    printf("mdaconvert input.ncs output.mda\n");
    printf("mdaconvert input.nrd output.mda --num_channels=32\n");
    printf("mdaconvert extract_time_chunk input.mda output.mda --t1=0 --t2=1e6\n");
    printf("mdaconvert extract_channels input.mda output.mda --channels=5,6,7,8,16-20\n");
}

#define MDAIO_MAX_DIMS 50
#define MDAIO_MAX_SIZE 1e15
#define MDAIO_TYPE_COMPLEX -1
#define MDAIO_TYPE_BYTE -2
#define MDAIO_TYPE_FLOAT32 -3
#define MDAIO_TYPE_INT16 -4
#define MDAIO_TYPE_INT32 -5
#define MDAIO_TYPE_UINT16 -6
#define MDAIO_TYPE_FLOAT64 -7
#define MDAIO_TYPE_UINT32 -8

QString mdaio_type_string(int dtype)
{
    if (dtype == MDAIO_TYPE_COMPLEX)
        return "complex";
    else if (dtype == MDAIO_TYPE_FLOAT32)
        return "float32";
    else if (dtype == MDAIO_TYPE_FLOAT64)
        return "float64";
    else if (dtype == MDAIO_TYPE_BYTE)
        return "byte";
    else if (dtype == MDAIO_TYPE_UINT16)
        return "uint16";
    else if (dtype == MDAIO_TYPE_UINT32)
        return "uint32";
    else if (dtype == MDAIO_TYPE_INT16)
        return "int16";
    else if (dtype == MDAIO_TYPE_INT32)
        return "int32";
    else
        return QString("unknown type %1").arg(dtype);
}

QString get_info_string(const DiskReadMda& X)
{
    QString str = QString("%1x%2").arg(X.N1()).arg(X.N2());
    if (X.N6() > 1) {
        str += QString("x%1x%2x%3x%4").arg(X.N3()).arg(X.N4()).arg(X.N5()).arg(X.N6());
    }
    else if (X.N5() > 1) {
        str += QString("x%1x%2x%3").arg(X.N3()).arg(X.N4()).arg(X.N5());
    }
    else if (X.N4() > 1) {
        str += QString("x%1x%2").arg(X.N3()).arg(X.N4());
    }
    else if (X.N3() > 1) {
        str += QString("x%1").arg(X.N3());
    }
    MDAIO_HEADER header = X.mdaioHeader();
    str += QString(" (%1)").arg(mdaio_type_string(header.data_type));
    return str;
}

QString get_data_type_string(int data_type)
{
    if (data_type == MDAIO_TYPE_BYTE)
        return "byte";
    if (data_type == MDAIO_TYPE_FLOAT32)
        return "float32";
    if (data_type == MDAIO_TYPE_FLOAT64)
        return "float64";
    if (data_type == MDAIO_TYPE_INT16)
        return "int16";
    if (data_type == MDAIO_TYPE_INT32)
        return "int32";
    if (data_type == MDAIO_TYPE_UINT16)
        return "uint16";
    if (data_type == MDAIO_TYPE_UINT32)
        return "uint32";
    return "";
}

QString get_json_header_string(const DiskReadMda& X)
{
    MDAIO_HEADER H = X.mdaioHeader();
    QJsonObject ret;
    ret["num_dims"] = H.num_dims;
    QJsonArray dims0;
    for (int i = 0; i < H.num_dims; i++) {
        dims0.push_back((long long)H.dims[i]);
    }
    ret["dims"] = dims0;
    ret["num_bytes_per_entry"] = H.num_bytes_per_entry;
    ret["header_size"] = (int)H.header_size;
    ret["data_type"] = H.data_type;
    ret["data_type_string"] = get_data_type_string(H.data_type);
    return QJsonDocument(ret).toJson(QJsonDocument::Indented);
}

bool extract_channels(QString input_fname, QString output_fname, const QMap<QString, QVariant>& params)
{
    DiskReadMda X(input_fname);
    Mda Y;

    QStringList channels_str = params["channels"].toString().split(",", QString::SkipEmptyParts);
    QList<int> channels = MLUtil::stringListToIntList(channels_str);
    int M2 = channels.count();
    Y.allocate(M2, X.N2());
    for (bigint t = 0; t < X.N2(); t++) {
        for (int m = 0; m < M2; m++) {
            Y.set(X.value(channels[m] - 1, t), m, t);
        }
    }

    MDAIO_HEADER H = X.mdaioHeader();
    if (H.data_type == MDAIO_TYPE_FLOAT32)
        return Y.write32(output_fname);
    else if (H.data_type == MDAIO_TYPE_FLOAT64)
        return Y.write64(output_fname);
    else if (H.data_type == MDAIO_TYPE_INT16)
        return Y.write16i(output_fname);
    else if (H.data_type == MDAIO_TYPE_UINT16)
        return Y.write16ui(output_fname);
    else if (H.data_type == MDAIO_TYPE_INT32)
        return Y.write16i(output_fname);
    else if (H.data_type == MDAIO_TYPE_UINT32)
        return Y.write16ui(output_fname);
    else if (H.data_type == MDAIO_TYPE_BYTE)
        return Y.write8(output_fname);
    else {
        qWarning() << "Unexpected data type: " << H.data_type;
        return false;
    }
}

bool extract_time_chunk(QString input_fname, QString output_fname, const QMap<QString, QVariant>& params)
{
    DiskReadMda X(input_fname);
    Mda Y;

    bigint t1 = params["t1"].toDouble(); //to double to handle scientific notation
    bigint t2 = params["t2"].toDouble();

    if (!X.readChunk(Y, 0, t1, X.N1(), t2 - t1 + 1)) {
        qWarning() << "Problem reading chunk.";
        return false;
    }
    MDAIO_HEADER H = X.mdaioHeader();
    if (H.data_type == MDAIO_TYPE_FLOAT32)
        return Y.write32(output_fname);
    else if (H.data_type == MDAIO_TYPE_FLOAT64)
        return Y.write64(output_fname);
    else if (H.data_type == MDAIO_TYPE_INT16)
        return Y.write16i(output_fname);
    else if (H.data_type == MDAIO_TYPE_UINT16)
        return Y.write16ui(output_fname);
    else if (H.data_type == MDAIO_TYPE_INT32)
        return Y.write16i(output_fname);
    else if (H.data_type == MDAIO_TYPE_UINT32)
        return Y.write16ui(output_fname);
    else if (H.data_type == MDAIO_TYPE_BYTE)
        return Y.write8(output_fname);
    else {
        qWarning() << "Unexpected data type: " << H.data_type;
        return false;
    }
}
