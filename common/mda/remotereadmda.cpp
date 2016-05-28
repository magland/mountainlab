/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/2/2016
*******************************************************/

#include "remotereadmda.h"
#include "diskreadmda.h"
#include <QStringList>
#include <QDir>
#include <QDateTime>
#include <taskprogress.h>
#include "cachemanager.h"
#include "msmisc.h"

#define REMOTE_READ_MDA_CHUNK_SIZE 5e5

struct RemoteReadMdaInfo {
    long N1, N2, N3;
    QString checksum;
    QDateTime file_last_modified;
};

class RemoteReadMdaPrivate {
public:
    RemoteReadMda* q;

    QString m_path;
    RemoteReadMdaInfo m_info;
    bool m_info_downloaded;
    HaltAgent* m_halt_agent;

    void download_info_if_needed();
    QString download_chunk_at_index(long ii);
};

RemoteReadMda::RemoteReadMda(const QString& path)
{
    d = new RemoteReadMdaPrivate;
    d->q = this;
    d->m_halt_agent = 0;
    this->setPath(path);
}

RemoteReadMda::RemoteReadMda(const RemoteReadMda& other)
{
    d = new RemoteReadMdaPrivate;
    d->q = this;
    d->m_halt_agent = 0;
    this->setPath(other.d->m_path);
}

void RemoteReadMda::operator=(const RemoteReadMda& other)
{
    /// TODO should I copy the halt_agent?
    this->setPath(other.d->m_path);
}

RemoteReadMda::~RemoteReadMda()
{
    delete d;
}

void RemoteReadMda::setPath(const QString& path)
{
    d->m_info.N1 = d->m_info.N2 = d->m_info.N3 = 0;
    d->m_info_downloaded = false;
    d->m_path = path;
}

QString RemoteReadMda::path() const
{
    return d->m_path;
}

void RemoteReadMda::setHaltAgent(HaltAgent* halt_agent)
{
    d->m_halt_agent = halt_agent;
}

long RemoteReadMda::N1()
{
    d->download_info_if_needed();
    return d->m_info.N1;
}

long RemoteReadMda::N2()
{
    d->download_info_if_needed();
    return d->m_info.N2;
}

long RemoteReadMda::N3()
{
    d->download_info_if_needed();
    return d->m_info.N3;
}

QDateTime RemoteReadMda::fileLastModified()
{
    d->download_info_if_needed();
    return d->m_info.file_last_modified;
}

bool RemoteReadMda::readChunk(Mda& X, long i, long size) const
{
    /// TODO handle both size=8 and 4 bytes
    int datatype_size=8;
    double size_mb=size*datatype_size*1.0/1e6;
    TaskProgress task;
    if (size_mb>0.5) {
        task.setLabel(QString("Downloading array chunk: %1 MB").arg(size_mb));
    }
    //read a chunk of the remote array considered as a 1D array
    X.allocate(size, 1); //allocate the output array
    double* Xptr = X.dataPtr(); //pointer to the output data
    long ii1 = i; //start index of the remote array
    long ii2 = i + size - 1; //end index of the remote array
    long jj1 = ii1 / REMOTE_READ_MDA_CHUNK_SIZE; //start chunk index of the remote array
    long jj2 = ii2 / REMOTE_READ_MDA_CHUNK_SIZE; //end chunk index of the remote array
    if (jj1 == jj2) { //in this case there is only one chunk we need to worry about
        task.setProgress(0.5);
        QString fname = d->download_chunk_at_index(jj1); //download the single chunk
        if (fname.isEmpty())
            return false;
        DiskReadMda A(fname);
        A.readChunk(X, ii1 - jj1 * REMOTE_READ_MDA_CHUNK_SIZE, size); //starting reading at the offset of ii1 relative to the start index of the chunk
        return true;
    } else {
        for (long jj = jj1; jj <= jj2; jj++) { //otherwise we need to step through the chunks
            task.setProgress((jj-jj1+0.5)/(jj2-jj1+1));
            if ((d->m_halt_agent) && (d->m_halt_agent->stopRequested())) {
                //X = Mda(); //maybe it's better to return the right size.
                return false;
            }
            QString fname = d->download_chunk_at_index(jj); //download the chunk at index jj
            if (fname.isEmpty())
                return false;
            DiskReadMda A(fname);
            if (jj == jj1) { //case 1/3, this is the first chunk
                Mda tmp;
                long size0 = (jj1 + 1) * REMOTE_READ_MDA_CHUNK_SIZE - ii1; //the size is going to be the difference between ii1 and the start index of the next chunk
                A.readChunk(tmp, ii1 - jj1 * REMOTE_READ_MDA_CHUNK_SIZE, size0); //again we start reading at the offset of ii1 relative to the start index of the chunk
                double* tmp_ptr = tmp.dataPtr(); //copy the data directly from tmp_ptr to Xptr
                long b = 0;
                for (long a = 0; a < size0; a++) {
                    Xptr[b] = tmp_ptr[a];
                    b++;
                }
            } else if (jj == jj2) { //case 2/3, this is the last chunk
                Mda tmp;
                long size0 = ii2 + 1 - jj2 * REMOTE_READ_MDA_CHUNK_SIZE; //the size is going to be the difference between the start index of the last chunk and ii2+1
                A.readChunk(tmp, 0, size0); //we start reading at position zero
                double* tmp_ptr = tmp.dataPtr();
                //copy the data to the last part of X
                long b = size - size0;
                for (long a = 0; a < size0; a++) {
                    Xptr[b] = tmp_ptr[a];
                    b++;
                }
            }

            /*
             ///this was the old code, which was wrong!!!!!!!!! -- fixed on 5/16/2016
            else if (jj == jj2) {
                Mda tmp;
                long size0 = ii2 + 1 - jj2 * REMOTE_READ_MDA_CHUNK_SIZE;
                A.readChunk(tmp, REMOTE_READ_MDA_CHUNK_SIZE - size0, size0);
                double* tmp_ptr = tmp.dataPtr();
                long b = jj2 * REMOTE_READ_MDA_CHUNK_SIZE - ii1;
                for (long a = 0; a < size0; a++) {
                    Xptr[b] = tmp_ptr[a];
                    b++;
                }
            }
            */
            else { //case 3/3, this is a middle chunk
                Mda tmp;
                A.readChunk(tmp, 0, REMOTE_READ_MDA_CHUNK_SIZE); //read the entire chunk, because we'll use it all
                double* tmp_ptr = tmp.dataPtr();
                long b = jj * REMOTE_READ_MDA_CHUNK_SIZE - ii1; //we start writing at the offset between the start index of the chunk and the start index
                for (long a = 0; a < REMOTE_READ_MDA_CHUNK_SIZE; a++) {
                    Xptr[b] = tmp_ptr[a];
                    b++;
                }
            }
        }
        return true;
    }
}

#ifdef USE_REMOTE_MDA
void RemoteReadMdaPrivate::download_info_if_needed()
{
    if (m_info_downloaded)
        return;
    m_info_downloaded = true;
    //QString url=file_url_for_remote_path(m_path);
    QString url = m_path;
    QString url2 = url + "?a=info&output=text";
    QString txt = http_get_text(url2);
    QStringList lines = txt.split("\n");
    QStringList sizes = lines.value(0).split(",");
    m_info.N1 = sizes.value(0).toLong();
    m_info.N2 = sizes.value(1).toLong();
    m_info.N3 = sizes.value(2).toLong();
    m_info.checksum = lines.value(1);
    m_info.file_last_modified = QDateTime::fromMSecsSinceEpoch(lines.value(2).toLong());
}
#endif

/// TODO only download float64 when necessary!
QString RemoteReadMdaPrivate::download_chunk_at_index(long ii)
{
    long Ntot = m_info.N1 * m_info.N2 * m_info.N3;
    long size = REMOTE_READ_MDA_CHUNK_SIZE;
    if (ii * REMOTE_READ_MDA_CHUNK_SIZE + size > Ntot) {
        size = Ntot - ii * REMOTE_READ_MDA_CHUNK_SIZE;
    }
    if (size <= 0)
        return "";
    if (m_info.checksum.isEmpty())
        return "";
    QString file_name = m_info.checksum + "-" + QString("%1-%2").arg(REMOTE_READ_MDA_CHUNK_SIZE).arg(ii);
    QString fname = CacheManager::globalInstance()->makeLocalFile(file_name, CacheManager::ShortTerm);
    if (QFile::exists(fname))
        return fname;
    QString url = m_path;
    QString url0 = url + QString("?a=readChunk&output=text&index=%1&size=%2&datatype=float64").arg((long)(ii * REMOTE_READ_MDA_CHUNK_SIZE)).arg(size);
    QString binary_url = http_get_text(url0).trimmed();
    if (binary_url.isEmpty())
        return "";

    //the following is ugly
    int ind = m_path.indexOf("/mdaserver");
    if (ind > 0) {
        binary_url = m_path.mid(0, ind) + "/mdaserver/" + binary_url;
    }

    QString mda_fname = http_get_binary_file(binary_url);
    if (mda_fname.isEmpty())
        return "";
    DiskReadMda tmp(mda_fname);
    if (tmp.totalSize() != size) {
        qWarning() << "Unexpected total size problem: " << tmp.totalSize() << size;
        QFile::remove(mda_fname);
        return "";
    }
    QFile::rename(mda_fname, fname);
    return fname;
}

void unit_test_remote_read_mda()
{
    QString url = "http://localhost:8000/firings.mda";
    RemoteReadMda X(url);
    qDebug() << X.N1() << X.N2() << X.N3();
    Mda chunk;
    X.readChunk(chunk, 0, 100);
    for (int j = 0; j < 10; j++) {
        qDebug() << j << chunk.value(j);
    }
}

#include "diskreadmda.h"
void unit_test_remote_read_mda_2(const QString& path)
{
    // run this test by calling
    // > mountainview unit_test remotereadmda2

    DiskReadMda X(path);
    for (int j = 0; j < 20; j++) {
        printf("%d: ", j);
        for (int i = 0; i < X.N1(); i++) {
            printf("%g, ", X.value(i, j));
        }
        printf("\n");
    }
}
