#include "msmisc.h"
#include <math.h>
#include <QDateTime>
#include <QDir>
#include <QCryptographicHash>
#include <QThread>
#include "cachemanager.h"

//I have to temporarily put this code in to get code completion to work in QtCreator
//even though DEFINES+=USE_NETWORK is in the .pro file.
//#define USE_NETWORK

#ifdef USE_NETWORK
#include <QEventLoop>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QTemporaryFile>
#include "taskprogress.h"
#endif

#include "textfile.h"

double compute_min(const QList<double>& X)
{
    double ret = X.value(0);
    for (int i = 0; i < X.count(); i++)
        if (X[i] < ret)
            ret = X[i];
    return ret;
}

double compute_max(const QList<double>& X)
{
    double ret = X.value(0);
    for (int i = 0; i < X.count(); i++)
        if (X[i] > ret)
            ret = X[i];
    return ret;
}

int compute_max(const QList<int>& X)
{
    int ret = X.value(0);
    for (int i = 0; i < X.count(); i++)
        if (X[i] > ret)
            ret = X[i];
    return ret;
}

long compute_max(const QList<long>& X)
{
    long ret = X.value(0);
    for (int i = 0; i < X.count(); i++)
        if (X[i] > ret)
            ret = X[i];
    return ret;
}

double compute_norm(long N, double* X)
{
    double sumsqr = 0;
    for (long i = 0; i < N; i++)
        sumsqr += X[i] * X[i];
    return sqrt(sumsqr);
}

Mda compute_mean_clip(Mda& clips)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    Mda ret;
    ret.allocate(M, T);
    int aaa = 0;
    for (int i = 0; i < L; i++) {
        int bbb = 0;
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.set(ret.get(bbb) + clips.get(aaa), bbb);
                aaa++;
                bbb++;
            }
        }
    }
    if (L) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.set(ret.get(m, t) / L, m, t);
            }
        }
    }
    return ret;
}

double compute_mean(const QList<double>& X)
{
    double sum = 0;
    for (int i = 0; i < X.count(); i++)
        sum += X[i];
    if (X.count())
        sum /= X.count();
    return sum;
}

double compute_stdev(const QList<double>& X)
{
    double sumsqr = 0;
    for (int i = 0; i < X.count(); i++)
        sumsqr += X[i] * X[i];
    double sum = 0;
    for (int i = 0; i < X.count(); i++)
        sum += X[i];
    int ct = X.count();
    if (ct >= 2) {
        return sqrt((sumsqr - sum * sum / ct) / (ct - 1));
    } else
        return 0;
}
Mda grab_clips_subset(Mda& clips, const QList<int>& inds)
{
    int M = clips.N1();
    int T = clips.N2();
    int LLL = inds.count();
    Mda ret;
    ret.allocate(M, T, LLL);
    for (int i = 0; i < LLL; i++) {
        long aaa = i * M * T;
        long bbb = inds[i] * M * T;
        for (int k = 0; k < M * T; k++) {
            ret.set(clips.get(bbb), aaa);
            aaa++;
            bbb++;
        }
    }
    return ret;
}

double compute_max(long N, double* X)
{
    if (N == 0)
        return 0;
    double ret = X[0];
    for (long i = 0; i < N; i++) {
        if (X[i] > ret)
            ret = X[i];
    }
    return ret;
}

QString get_temp_fname()
{
    return CacheManager::globalInstance()->makeLocalFile();
    //long rand_num = qrand() + QDateTime::currentDateTime().toMSecsSinceEpoch();
    //return QString("%1/MdaClient_%2.tmp").arg(QDir::tempPath()).arg(rand_num);
}

QString abbreviate(const QString& str, int len1, int len2)
{
    if (str.count() <= len1 + len2 + 20)
        return str;
    return str.mid(0, len1) + "...\n...\n..." + str.mid(str.count() - len2);
}

#ifdef USE_NETWORK
QString http_get_binary_file(const QString& url)
{
    QTime timer;
    timer.start();
    QString fname = get_temp_fname();
    QNetworkAccessManager manager; // better make it a singleton
    QNetworkReply* reply = manager.get(QNetworkRequest(QUrl(url)));
    QEventLoop loop;
    QFile temp(fname);
    long num_bytes = 0;
    temp.open(QIODevice::WriteOnly);
    QObject::connect(reply, &QNetworkReply::readyRead, [&]() {
        QByteArray X=reply->readAll();
        temp.write(X);
        num_bytes+=X.count();
    });
    QObject::connect(reply, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();
    printf("RECEIVED BINARY (%d ms, %ld bytes) from %s\n", timer.elapsed(), num_bytes, url.toLatin1().data());

    TaskProgressAgent::globalInstance()->incrementQuantity("bytes_downloaded", num_bytes);

    return fname;
}

QString http_get_text(const QString& url)
{
    QTime timer;
    timer.start();
    QString fname = get_temp_fname();
    QNetworkAccessManager manager; // better make it a singleton
    QNetworkReply* reply = manager.get(QNetworkRequest(QUrl(url)));
    QEventLoop loop;
    QString ret;
    QObject::connect(reply, &QNetworkReply::readyRead, [&]() {
        ret+=reply->readAll();
    });
    QObject::connect(reply, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();
    printf("RECEIVED TEXT (%d ms, %d bytes) from GET %s\n", timer.elapsed(), ret.count(), url.toLatin1().data());
    QString str = abbreviate(ret, 200, 200);
    printf("%s\n", (str.toLatin1().data()));

    TaskProgressAgent::globalInstance()->incrementQuantity("bytes_downloaded", ret.count());

    return ret;
}
#else
QString http_get_binary_file(const QString& url)
{
    QString tmp_fname = get_temp_fname();
    QString cmd = QString("curl \"%1\" > %2").arg(url).arg(tmp_fname);
    qDebug() << cmd;
    int exit_code = system(cmd.toLatin1().data());
    if (exit_code != 0) {
        qWarning() << "Problem with system call: " + cmd;
        QFile::remove(tmp_fname);
        return "";
    }
    return tmp_fname;
}
QString http_get_text(const QString& url)
{
    QString tmp_fname = get_temp_fname();
    QString cmd = QString("curl \"%1\" > %2").arg(url).arg(tmp_fname);
    qDebug() << cmd;
    int exit_code = system(cmd.toLatin1().data());
    if (exit_code != 0) {
        qWarning() << "Problem with system call: " + cmd;
        QFile::remove(tmp_fname);
        return "";
    }
    QString ret = read_text_file(tmp_fname);
    QFile::remove(tmp_fname);
    qDebug() << "RESPONSE: " << ret;
    return ret;
}
#endif

QString compute_hash(const QString& str)
{
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(str.toLatin1());
    return QString(hash.result().toHex());
}
