#include "toolbuttonmenu.h"
#include <QEvent>
#include <QGraphicsOpacityEffect>
#include <QPropertyAnimation>
#include <QMenu>
#include <QTimer>
#include <QtDebug>

class FadingButton : public QToolButton {

public:

    const int interval = 2000;

    FadingButton(QWidget *parent = 0)
        : QToolButton(parent), m_state(Hidden) {
        m_timer.setInterval(interval);
        m_timer.setSingleShot(true);
        QGraphicsOpacityEffect* eff = new QGraphicsOpacityEffect(this);
        setGraphicsEffect(eff);
        eff->setOpacity(0);
        connect(&m_timer, &QTimer::timeout, [this]() {
            if (!isHovering()) {
                fadeOut();
            }
        });
        setAttribute(Qt::WA_Hover);
        m_anim = new QPropertyAnimation(eff, "opacity");
        m_anim->setDuration(250);
        connect(m_anim, &QAbstractAnimation::finished, [this](){
            if (m_state == FadingIn) {
                m_state = Shown;
            } else if (m_state == FadingOut) {
                m_state = Hidden;
            }
        });
        installEventFilter(this);
    }
    void fadeIn() {
        if (m_state == Hidden || m_state == FadingOut) {
            m_anim->setStartValue(0);
            m_anim->setEndValue(1);
            m_anim->start();
            m_state = FadingIn;
        }
        start();
    }
    void fadeOut() {
        if (m_state == Shown || m_state == FadingIn) {
            m_anim->setStartValue(1);
            m_anim->setEndValue(0);
            m_anim->start();
            m_state = FadingOut;
        }
        stop();
    }
    void forceStart() {
        m_hovering = false;
        start();
    }

    void start() {
        m_timer.start();

    }
    void stop() {
        m_timer.stop();
    }
    bool isRunning() const { return m_timer.isActive(); }
    bool isHovering() const { return m_hovering; }
protected:
    bool eventFilter(QObject *o, QEvent *e) {
        if (o != this) return QToolButton::eventFilter(o, e);

        if (e->type() == QEvent::HoverEnter || e->type() == QEvent::HoverMove) {
            m_hovering = true;
//            stop();
        } else if (e->type() == QEvent::HoverLeave) {
            m_hovering = false;
            start();
        }
        return QToolButton::eventFilter(o, e);
    }

private:
    QTimer m_timer;
    QPropertyAnimation *m_anim;
    bool m_hovering = false;
    enum { Hidden, FadingIn, Shown, FadingOut } m_state;

};

ToolButtonMenu::ToolButtonMenu(QObject* parent)
    : QObject(parent)
    , m_offset(QSize(4, 0))
{
}

void ToolButtonMenu::setOffset(const QSize& off)
{
    if (offset() == off)
        return;
    m_offset = off;
    emit offsetChanged(m_offset);
}

QToolButton* ToolButtonMenu::activateOn(QWidget* w)
{
    w->installEventFilter(this);
    w->setAttribute(Qt::WA_Hover);
    FadingButton* tb = new FadingButton(w);

    tb->setAutoRaise(true);
    connect(tb, SIGNAL(clicked()), this, SLOT(openMenu()));

    m_widgets.insert(w, tb);
    tb->setVisible(true);
    updateButtonPosition(w, tb);
    return tb;
}

bool ToolButtonMenu::eventFilter(QObject* o, QEvent* e)
{
    if (o->isWidgetType()) {
        QWidget* w = static_cast<QWidget*>(o);
        if (e->type() == QEvent::HoverEnter) {
            FadingButton* tb = toolButton(w);
            tb->fadeIn();
            tb->start();
        }
        else if (e->type() == QEvent::HoverEnter || e->type() == QEvent::HoverMove) {
            FadingButton* tb = toolButton(w);
            tb->fadeIn();
        }
        else if (e->type() == QEvent::HoverLeave) {
            FadingButton* tb = toolButton(w);
            tb->fadeOut();
        }
        else if (e->type() == QEvent::Resize) {
            updateButtonPosition(w);
        }
    }
    return false;
}

void ToolButtonMenu::updateButtonPosition(QWidget* w, QToolButton* tb)
{
    if (!tb)
        tb = toolButton(w);
    if (!w || !tb)
        return;
    QSize sh = tb->sizeHint();
    tb->setGeometry(w->width() - sh.width() - offset().width(), offset().height(),
                    sh.width(), sh.height());
}

FadingButton *ToolButtonMenu::toolButton(QWidget *w) const
{
    return m_widgets.value(w, Q_NULLPTR);
}


QWidget *ToolButtonMenu::widget(FadingButton *tb) const
{
    return m_widgets.key(tb, Q_NULLPTR);
}

void ToolButtonMenu::openMenu()
{
    FadingButton* tb = static_cast<FadingButton*>(qobject_cast<QToolButton*>(sender()));
    if (!tb)
        return;
    QWidget* w = widget(tb);
    if (!w)
        return;
    QMenu m;
    m.addActions(w->actions());
    m.exec(tb->mapToGlobal(tb->rect().bottomLeft()));
    if (tb->rect().contains(tb->mapFromGlobal(QCursor::pos())))
        tb->start();
    else
        tb->forceStart();
}
