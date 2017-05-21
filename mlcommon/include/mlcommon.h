/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/6/2016
*******************************************************/

#ifndef MLCOMMON_H
#define MLCOMMON_H

#include <QTextCodec>
#include <QDebug>
#include <QJsonValue>
#include <QByteArray>
#include <QJsonObject>

#include <math.h>

typedef int64_t bigint;

namespace TextFile {
QString read(const QString& fname, QTextCodec* codec = 0);
bool write(const QString& fname, const QString& txt, QTextCodec* codec = 0);
bool write_single_try(const QString& fname, const QString& txt, QTextCodec* codec = 0);
};

namespace MLUtil {
QString makeRandomId(int numchars = 10);
bool threadInterruptRequested();
bool inGuiThread();
QString tempPath();
QString mountainlabBasePath();
QString mlLogPath();
QString resolvePath(const QString& basepath, const QString& path);
void mkdirIfNeeded(const QString& path);
QString computeSha1SumOfFile(const QString& path);
QString computeSha1SumOfFileHead(const QString& path, bigint num_bytes);
QString computeSha1SumOfString(const QString& str);
QString computeSha1SumOfDirectory(const QString& path);
bool matchesFastChecksum(QString path, QString fcs);
QList<int> stringListToIntList(const QStringList& list);
QList<bigint> stringListToBigIntList(const QStringList& list);
QStringList intListToStringList(const QList<int>& list);
QJsonValue toJsonValue(const QByteArray& X);
QJsonValue toJsonValue(const QList<int>& X);
QJsonValue toJsonValue(const QVector<int>& X);
QJsonValue toJsonValue(const QVector<double>& X);
void fromJsonValue(QByteArray& X, const QJsonValue& val);
void fromJsonValue(QList<int>& X, const QJsonValue& val);
void fromJsonValue(QVector<int>& X, const QJsonValue& val);
void fromJsonValue(QVector<double>& X, const QJsonValue& val);
QByteArray readByteArray(const QString& path);
bool writeByteArray(const QString& path, const QByteArray& X);
QJsonObject mountainlabConfig();
QJsonValue configValue(const QString& group, const QString& key);
QString configResolvedPath(const QString& group, const QString& key);
QStringList configResolvedPathList(const QString& group, const QString& key);
QStringList toStringList(const QVariant& val); //val is either a string or a QVariantList
QJsonObject createPrvObject(const QString& file_or_dir_path);
QString locatePrv(const QJsonObject& obj, const QStringList& local_search_paths);
};

namespace MLCompute {
double min(const QVector<double>& X);
double max(const QVector<double>& X);
double sum(const QVector<double>& X);
double mean(const QVector<double>& X);
double median(const QVector<double>& X);
double stdev(const QVector<double>& X);
double norm(const QVector<double>& X);
double dotProduct(const QVector<double>& X1, const QVector<double>& X2);
double correlation(const QVector<double>& X1, const QVector<double>& X2);

double norm(const QVector<float>& X);
double dotProduct(const QVector<float>& X1, const QVector<float>& X2);

template <typename T>
T max(const QVector<T>& X);
double min(bigint N, const double* X);
double max(bigint N, const double* X);
double sum(bigint N, const double* X);
double mean(bigint N, const double* X);
double dotProduct(bigint N, const double* X1, const double* X2);
double norm(bigint N, const double* X);
double min(bigint N, const float* X);
double max(bigint N, const float* X);
double sum(bigint N, const float* X);
double mean(bigint N, const float* X);
double stdev(bigint N, const float* X);
double dotProduct(bigint N, const float* X1, const float* X2);
double correlation(bigint N, const float* X1, const float* X2);
double norm(bigint N, const float* X);
}

class CLParams {
public:
    CLParams(int argc, char* argv[]);
    QMap<QString, QVariant> named_parameters;
    QList<QString> unnamed_parameters;
    bool success;
    QString error_message;
};

template <typename T>
T MLCompute::max(const QVector<T>& X)
{
    return *std::max_element(X.constBegin(), X.constEnd());
}

QString locate_prv(const QJsonObject& obj);

/*
QString resolve_prv_object(const QJsonObject& obj, bool allow_downloads, bool allow_processing);
QString resolve_prv_file(const QString& prv_fname, bool allow_downloads, bool allow_processing);
bool resolve_prv_files(QMap<QString, QVariant>& command_line_params, bool allow_downloads, bool allow_processing);
bool prepare_prv_files(QMap<QString, QVariant>& command_line_params, bool allow_downloads, bool allow_processing);
*/

#endif // TEXTFILE_H
