/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 8/25/2016
*******************************************************/

#include "clustercontextmenuplugin.h"
#include <clustercontextmenuhandler.h>
#include <clusterpaircontextmenuhandler.h>

class ClusterContextMenuPluginPrivate {
public:
    ClusterContextMenuPlugin* q;
};

ClusterContextMenuPlugin::ClusterContextMenuPlugin()
{
    d = new ClusterContextMenuPluginPrivate;
    d->q = this;
}

ClusterContextMenuPlugin::~ClusterContextMenuPlugin()
{
    delete d;
}

QString ClusterContextMenuPlugin::name()
{
    return "ClusterContextMenu";
}

QString ClusterContextMenuPlugin::description()
{
    return "Context menu handlers for selected clusters or selected cluster pairs.";
}

void ClusterContextMenuPlugin::initialize(MVMainWindow* mw)
{
    mw->registerContextMenuHandler(new MVClusterContextMenuHandler(mw));
    //mw->registerContextMenuHandler(new MVClusterPairContextMenuHandler(mw));
    MVAbstractPlugin::initialize(mw);
}
