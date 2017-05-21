#ifndef MCVIEWFACTORIES_H
#define MCVIEWFACTORIES_H

#include "mccontext.h"

#include <mvabstractviewfactory.h>

/////////////////////////////////////////////////////////////////////////////
class ClusterDetail1Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    ClusterDetail1Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
private slots:
    //void openClipsForTemplate();
};

class ClusterDetail2Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    ClusterDetail2Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
private slots:
    //void openClipsForTemplate();
};

/////////////////////////////////////////////////////////////////////////////
class MVTimeSeriesView1Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVTimeSeriesView1Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
};

class MVTimeSeriesView2Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVTimeSeriesView2Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
};

class MVTimeSeriesViewIntersectFactory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVTimeSeriesViewIntersectFactory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
};

/////////////////////////////////////////////////////////////////////////////
class MVSpikeSpray1Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVSpikeSpray1Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
    //QList<QAction*> actions(const QMimeData& md) Q_DECL_OVERRIDE;
};

class MVSpikeSpray2Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVSpikeSpray2Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
    //QList<QAction*> actions(const QMimeData& md) Q_DECL_OVERRIDE;
};

/////////////////////////////////////////////////////////////////////////////
class MVPCAFeatures1Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVPCAFeatures1Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
    bool isEnabled(MVAbstractContext* context) const Q_DECL_OVERRIDE;
};

class MVPCAFeatures2Factory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    MVPCAFeatures2Factory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
    bool isEnabled(MVAbstractContext* context) const Q_DECL_OVERRIDE;
};

/////////////////////////////////////////////////////////////////////////////
class CompareClustersFactory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    CompareClustersFactory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
private slots:
    //void openClipsForTemplate();
};

#endif // MCVIEWFACTORIES_H
