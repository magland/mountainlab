/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef TABBER_H
#define TABBER_H

#include <QTabWidget>
#include <QString>

/**
 *  \class Tabber
 *  \brief Use this class to control two or more tab widgets so that the user
 * may move widgets from one to the other. For now user double clicks on a tab
 * to move it to the other, but in the future, drag/drop should be enabled. Also
 * we'd like user to be able to drag it off the widget to become undocked (floating).
 */

class TabberTabWidget;
class TabberPrivate;
class Tabber : public QObject {
    Q_OBJECT
public:
    friend class TabberPrivate;
    Tabber();
    virtual ~Tabber();
    ///Create a new tab widget that will be under the controll of this Tabber
    TabberTabWidget* createTabWidget(const QString& container_name);
    ///Add a widget to a particular container (e.g., tab widget) with a given tab label
    void addWidget(const QString& container_name, const QString& label, QWidget* W);
    ///Add a widget to a particular container with a given tab label
    void addWidget(TabberTabWidget* TW, const QString& label, QWidget* W);
    ///The name of the last clicked container, and it controls where the next addWidget will go
    QString currentContainerName();
    ///Set the current container by name (see currentContainerName)
    void setCurrentContainerName(const QString& container_name);
    ///Set the current container
    void setCurrentContainer(TabberTabWidget* TW);
    ///Move the current container to another one -- or if there are only to, the "other" one. Useful for popping up a view on a different view.
    void switchCurrentContainer();
    ///Retrieve all widgets that have been added and have not been removed and deleted
    QList<QWidget*> allWidgets();
private
slots:
    void slot_tab_close_requested(int index);
    void slot_tab_bar_clicked(int index);
    void slot_tab_bar_double_clicked(int index);
    void slot_widget_destroyed(QObject* obj);

private:
    TabberPrivate* d;
};

class TabberTabWidgetPrivate;
class TabberTabWidget : public QTabWidget {
    Q_OBJECT
public:
    friend class TabberTabWidgetPrivate;
    TabberTabWidget();
    virtual ~TabberTabWidget();

private:
    TabberTabWidgetPrivate* d;
};

#endif // TABBER_H
