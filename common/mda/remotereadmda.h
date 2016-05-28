/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/2/2016
*******************************************************/

#ifndef REMOTEREADMDA_H
#define REMOTEREADMDA_H

#include <QString>
#include "mda.h"
#include "haltagent.h"

class RemoteReadMdaPrivate;
/**
 * @brief Do not use this class directly -- it is used by DiskReadMda
 */
class RemoteReadMda {
public:
    friend class RemoteReadMdaPrivate;
    RemoteReadMda(const QString& path = "");
    RemoteReadMda(const RemoteReadMda& other);
    void operator=(const RemoteReadMda& other);
    virtual ~RemoteReadMda();

    void setPath(const QString& path);
    QString path() const;
    void setHaltAgent(HaltAgent* halt_agent);

    long N1();
    long N2();
    long N3();
    QDateTime fileLastModified();

    ///Retrieve a chunk of the vectorized data of size 1xN starting at position i
    bool readChunk(Mda& X, long i, long size) const;

private:
    RemoteReadMdaPrivate* d;
};

void unit_test_remote_read_mda_2(const QString& path);
void unit_test_remote_read_mda();

#endif // REMOTEREADMDA_H
