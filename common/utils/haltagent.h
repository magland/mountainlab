/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/26/2016
*******************************************************/

#ifndef HALTAGENT_H
#define HALTAGENT_H

#include <QObject>

class HaltAgentPrivate;
class HaltAgent
{
public:
    virtual bool stopRequested()=0;
};

#endif // HALTAGENT_H
