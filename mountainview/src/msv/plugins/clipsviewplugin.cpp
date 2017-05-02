/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 8/24/2016
*******************************************************/

#include "clipsviewplugin.h"

#include <QAction>
#include <mvclipswidget.h>
#include <mvcontext.h>

class ClipsViewPluginPrivate {
public:
    ClipsViewPlugin* q;
};

ClipsViewPlugin::ClipsViewPlugin()
{
    d = new ClipsViewPluginPrivate;
    d->q = this;
}

ClipsViewPlugin::~ClipsViewPlugin()
{
    delete d;
}

QString ClipsViewPlugin::name()
{
    return "ClipsView";
}

QString ClipsViewPlugin::description()
{
    return "";
}

void ClipsViewPlugin::initialize(MVMainWindow* mw)
{
    mw->registerViewFactory(new ClipsViewFactory(mw));
    MVAbstractPlugin::initialize(mw);
}

QList<QAction*> ClipsViewPlugin::actions(const QMimeData& md)
{
    (void)md;
    return QList<QAction*>();
    /*
    QList<QAction*> actions;
    if (!md.data("x-mv-main").isEmpty()) {
        QAction* action = new QAction("Clips", 0);
        MVMainWindow* mw = this->mainWindow();
        QObject::connect(action, &QAction::triggered, [mw]() {
            mw->openView("open-clips-view");
        });
        actions << action;
    }
    return actions;
    */
}

ClipsViewFactory::ClipsViewFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString ClipsViewFactory::id() const
{
    return QStringLiteral("open-clips-view");
}

QString ClipsViewFactory::name() const
{
    return tr("Clips");
}

QString ClipsViewFactory::title() const
{
    return tr("Clips");
}

MVAbstractView* ClipsViewFactory::createView(MVAbstractContext* context)
{
    MVContext* c = qobject_cast<MVContext*>(context);
    Q_ASSERT(c);
    MVClipsWidget* X = new MVClipsWidget(c);
    X->setLabelsToUse(c->selectedClusters());
    return X;
}
