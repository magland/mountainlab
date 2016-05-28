#include "testMdaIO.h"
#include "mdaio.h"
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <QTemporaryFile>
#include <QDebug>

class TempFile {
public:
    TempFile()
    {
        if (this->file.open()) {
            this->handle = fdopen(this->file.handle(), "r+b");
        } else {
            this->handle = 0;
        }
    }
    ~TempFile()
    {
        if (this->handle) {
            const QString path = this->file.fileName();
            this->file.close();
            fclose(this->handle);
            if (QFile::exists(path))
                QFile::remove(path);
        }
    }
    bool isOpen() const
    {
        return this->file.isOpen();
    }
    void rewind()
    {
        ::rewind(this->handle);
    }
    operator FILE*()
    {
        return this->handle;
    }

private:
    QTemporaryFile file;
    FILE* handle;
};

void TestMdaIO::testHeader()
{
    struct MDAIO_HEADER header, headerCopy;
    header.data_type = MDAIO_TYPE_FLOAT32;
    const size_t dimCount = sizeof(header.dims) / sizeof(header.dims[0]);
    for (size_t i = 0; i < dimCount; ++i)
        header.dims[i] = 1 + i % 5;
    header.num_dims = dimCount / 2;
    header.num_bytes_per_entry = 4;
    mda_copy_header(&headerCopy, &header);
    QVERIFY(std::memcmp(&header, &headerCopy, sizeof(header)) == 0);
    TempFile file;
    QVERIFY(file.isOpen());
    QVERIFY(mda_write_header(&headerCopy, file));
    std::memset(&header, 0, sizeof(header));
    file.rewind();
    QVERIFY(mda_read_header(&header, file));
    QCOMPARE(header.data_type, headerCopy.data_type);
    QCOMPARE(header.num_dims, headerCopy.num_dims);
    QCOMPARE(header.num_bytes_per_entry, headerCopy.num_bytes_per_entry);
    QCOMPARE(header.header_size, headerCopy.header_size);
    for (int i = 0; i < header.num_dims; ++i)
        QCOMPARE(header.dims[i], headerCopy.dims[i]);
}

template <typename T>
std::vector<T> gen_data(const long size)
{
    std::vector<T> result(size);
    for (long i = 0; i < size; ++i) {
        result[i] = static_cast<T>(qrand() - RAND_MAX / 2);
    }
    return result;
}

template <typename T>
class _mdaiofunc {
public:
    typedef long (*func_type)(T*, struct MDAIO_HEADER*, long, FILE*);
    _mdaiofunc(func_type f)
        : _f(f)
    {
    }
    long operator()(T* data, struct MDAIO_HEADER* header, long size, FILE* file)
    {
        return this->_f(data, header, size, file);
    }

private:
    func_type _f;
};

template <typename T>
void test_type(const int typeId, _mdaiofunc<T> readFunc, _mdaiofunc<T> writeFunc)
{
    struct MDAIO_HEADER header;
    std::memset(&header, 0, sizeof(header));
    header.data_type = typeId;
    const long minimumSize = 8;
    const long size = qrand() % 100 + minimumSize;
    std::vector<T> data = gen_data<T>(size);
    TempFile file;
    QVERIFY(file.isOpen());
    QVERIFY(writeFunc(&data[0], &header, size, file));
    std::vector<T> result(size);
    file.rewind();
    QVERIFY(readFunc(&result[0], &header, size, file));
    bool compareResult = true;
    if (data != result) {
        compareResult = false;
        qDebug() << "---> data written to file:";
        std::copy(data.begin(), data.end(), std::ostream_iterator<T>(std::cout, ", "));
        qDebug() << "\n---> data read from file:";
        std::copy(result.begin(), result.end(), std::ostream_iterator<T>(std::cout, ", "));
    }
    QVERIFY(compareResult);
}

void TestMdaIO::testByte()
{
    qsrand(time(0));
    test_type<unsigned char>(MDAIO_TYPE_BYTE, mda_read_byte, mda_write_byte);
}

void TestMdaIO::testFloat32()
{
    test_type<float>(MDAIO_TYPE_FLOAT32, mda_read_float32, mda_write_float32);
}

void TestMdaIO::testInt16()
{
    test_type<int16_t>(MDAIO_TYPE_INT16, mda_read_int16, mda_write_int16);
}

void TestMdaIO::testInt32()
{
    test_type<int32_t>(MDAIO_TYPE_INT32, mda_read_int32, mda_write_int32);
}

void TestMdaIO::testUInt16()
{
    test_type<uint16_t>(MDAIO_TYPE_UINT16, mda_read_uint16, mda_write_uint16);
}

void TestMdaIO::testFloat64()
{
    test_type<double>(MDAIO_TYPE_FLOAT64, mda_read_float64, mda_write_float64);
}
