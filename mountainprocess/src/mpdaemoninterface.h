/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/2/2016
*******************************************************/

#ifndef MPDAEMONINTERFACE_H
#define MPDAEMONINTERFACE_H

#include <QJsonObject>


class MPDaemonInterfacePrivate;
class MPDaemonInterface
{
public:
    friend class MPDaemonInterfacePrivate;
    MPDaemonInterface();
    virtual ~MPDaemonInterface();
    bool start();
    bool stop();
    QJsonObject getInfo();
private:
    MPDaemonInterfacePrivate *d;
};

#endif // MPDAEMONINTERFACE_H
