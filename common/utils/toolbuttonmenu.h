#ifndef TOOLBUTTONMENU_H
#define TOOLBUTTONMENU_H

#include <QToolButton>
#include <QHash>
class QTimer;

class FadingButton;
class ToolButtonMenu : public QObject {
    Q_OBJECT
    Q_PROPERTY(QSize offset READ offset WRITE setOffset NOTIFY offsetChanged)
public:
    ToolButtonMenu(QObject* parent = 0);
    inline QSize offset() const { return m_offset; }
    void setOffset(const QSize& off);

    QToolButton* activateOn(QWidget* w);
signals:
    void offsetChanged(QSize);

protected:
    bool eventFilter(QObject* o, QEvent* e);
    void updateButtonPosition(QWidget* w, QToolButton* tb = 0);
    FadingButton *toolButton(QWidget *w) const;
    QWidget *widget(FadingButton *tb) const;

private slots:
    void openMenu();

private:

    QHash<QWidget*, FadingButton *> m_widgets;
    QSize m_offset;
};

#endif // TOOLBUTTONMENU_H
