/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/30/2016
*******************************************************/

#include "taskprogressview.h"
#include "taskprogress.h"
#include <QTreeWidget>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QStyledItemDelegate>
#include <QPainter>
#include <QtDebug>
#include <QPlainTextEdit>
#include <QDialogButtonBox>
#include <QDialog>
#include <QShortcut>
#include <QApplication>
#include <QClipboard>
#include <QDesktopWidget>
#include <QStandardItemModel>
#include <QTimer>
#include <QAction>
#include <QMenu>
#include <QSortFilterProxyModel>
#include <QScrollBar>
#include <QToolBar>
#include <QToolButton>

class TaskProgressViewDelegate : public QStyledItemDelegate {
public:
    TaskProgressViewDelegate(QObject* parent = 0)
        : QStyledItemDelegate(parent)
    {
    }
    QSize sizeHint(const QStyleOptionViewItem& option, const QModelIndex& index) const
    {
        if (index.parent().isValid())
            return QStyledItemDelegate::sizeHint(option, index);
        QSize sh = QStyledItemDelegate::sizeHint(option, index);
        sh.setHeight(sh.height() * 2);
        return sh;
    }
    void paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index) const
    {
        //        if (index.internalId() != 0xDEADBEEF || index.column() != 0) {
        if (index.parent().isValid() || index.column() != 0) {
            QStyledItemDelegate::paint(painter, option, index);
            return;
        }
        QSet<QString> tags = index.data(TaskManager::TaskProgressModel::TagsSetRole).value<QSet<QString> >();
        QIcon icon;
        if (tags.contains("calculate")) {
            QPixmap px(":/images/calculator.png");
            icon.addPixmap(px);
        }
        else if (tags.contains("download")) {
            QPixmap px(":/images/down.png");
            icon.addPixmap(px);
        }
        else if (tags.contains("process")) {
            QPixmap px(":/images/process.png");
            icon.addPixmap(px);
        }
        QStyleOptionViewItem opt = option;
        opt.text = "";
        opt.displayAlignment = Qt::AlignTop | Qt::AlignLeft;
        opt.decorationPosition = QStyleOptionViewItem::Left;
        opt.decorationAlignment = Qt::AlignTop | Qt::AlignLeft;
        opt.features |= QStyleOptionViewItem::HasDecoration;
        opt.icon = icon;
        QStyledItemDelegate::paint(painter, opt, index);
        qreal progress = index.data(Qt::UserRole).toDouble();
        if (progress < 1.0) {
            QPen p = painter->pen();
            QFont f = painter->font();
            p.setColor((option.state & QStyle::State_Selected) ? Qt::white : Qt::darkGray);
            f.setPointSize(f.pointSize() - 3);
            QFontMetrics smallFm(f);
            int elapsedWidth = smallFm.width("MMMMM");
            QStyleOptionProgressBar progOpt;
            progOpt.initFrom(opt.widget);
            progOpt.rect = opt.rect;
            progOpt.minimum = 0;
            progOpt.maximum = 100;
            progOpt.progress = progress * 100;
            progOpt.rect.setTop(progOpt.rect.center().y());
            progOpt.rect.adjust(4 + elapsedWidth + 4, 2, -4, -2);
            if (option.widget) {
                option.widget->style()->drawControl(QStyle::CE_ProgressBar, &progOpt, painter, option.widget);
            }
            painter->save();
            painter->setPen(p);
            painter->setFont(f);
            QRect r = opt.rect;
            r.setTop(option.rect.center().y());
            r.adjust(4, 2, -4, -2);
            r.setRight(r.left() + elapsedWidth);

            qreal duration = index.data(Qt::UserRole + 1).toDateTime().msecsTo(QDateTime::currentDateTime()) / 1000.0;
            if (duration < 100)
                painter->drawText(r, Qt::AlignRight | Qt::AlignVCenter, QString("%1s").arg(duration, 0, 'f', 2));
            else
                painter->drawText(r, Qt::AlignRight | Qt::AlignVCenter, QString("%1s").arg(qRound(duration)));
            painter->restore();
        }
        else {
            painter->save();
            QPen p = painter->pen();
            QFont f = painter->font();
            p.setColor((option.state & QStyle::State_Selected) ? Qt::white : Qt::darkGray);
            f.setPointSize(f.pointSize() - 2);
            painter->setPen(p);
            painter->setFont(f);
            QRect r = option.rect;
            r.setTop(option.rect.center().y());
            r.adjust(4, 2, -4, -2);

            qreal duration = index.data(Qt::UserRole + 1).toDateTime().msecsTo(index.data(Qt::UserRole + 2).toDateTime()) / 1000.0;
            painter->drawText(r, Qt::AlignLeft | Qt::AlignVCenter, QString("Completed in %1s").arg(duration));
            painter->restore();
        }
    }
};

class TaskProgressViewModeProxy : public QSortFilterProxyModel {
public:
    enum CompleteTasksMode {
        CTM_Shown,
        CTM_Hidden,
        CTM_HiddenIfOlderThan
    };
    enum LogsMode {
        LM_Shown,
        LM_Hidden
    };

    TaskProgressViewModeProxy(QObject *parent = 0) 
        : QSortFilterProxyModel(parent) 
    {
        m_showLogs = new QAction(this);
        m_showLogs->setText("Display logs");
        m_showLogs->setCheckable(true);
        m_showLogs->setChecked(true);
        m_showLogs->setIcon(QIcon::fromTheme("view-history"));
        connect(m_showLogs, &QAction::toggled, [this](bool t) {
            setLogsMode(t ? LM_Shown : LM_Hidden);
        });

//        connect(m_showLogs, SIGNAL(toggled(bool)), this, SLOT(updateLogsMode(bool)));
    }

    void setCompleteTasksMode(CompleteTasksMode m)
    {
        if (m == m_cMode)
            return;
        m_cMode = m;
        invalidateFilter();
    }

    void setLogsMode(LogsMode m)
    {
        if (m == m_lMode)
            return;
        m_lMode = m;
        invalidateFilter();
    }

    void setCompleteTasksModeThreshold(unsigned int secs)
    {
        if (secs == m_ctmThreshold)
            return;
        m_ctmThreshold = secs;
        if (completeTasksMode() == CTM_HiddenIfOlderThan)
            invalidateFilter();
    }

    LogsMode logsMode() const { return m_lMode; }
    CompleteTasksMode completeTasksMode() const { return m_cMode; }
    unsigned int completeTasksModeThreshold() const { return m_ctmThreshold; }

    QAction* showLogsAction() const
    {
        return m_showLogs;
    }

protected:
    bool filterAcceptsRow(int source_row, const QModelIndex& source_parent) const
    {
        TaskManager::TaskProgressModel* mdl = taskModel();
        if (!mdl)
            return true;
        QModelIndex sourceIndex = mdl->index(source_row, 0, source_parent);
        if (!sourceIndex.isValid())
            return false;
        bool isTask = mdl->isTask(sourceIndex);
        if (!isTask && logsMode() == LM_Hidden)
            return false;
        if (isTask) {
            if (!mdl->isActive(sourceIndex)) {
                if (completeTasksMode() == CTM_Shown)
                    return true;
                if (completeTasksMode() == CTM_Hidden)
                    return false;
                if (completeTasksMode() == CTM_HiddenIfOlderThan
                    && !mdl->isCompletedWithin(sourceIndex, completeTasksModeThreshold()))
                    return false;
            }
        }
        return true;
    }

    TaskManager::TaskProgressModel* taskModel() const
    {
        return (TaskManager::TaskProgressModel*)sourceModel();
    }

private:
    CompleteTasksMode m_cMode = CTM_Shown;
    LogsMode m_lMode = LM_Shown;
    unsigned int m_ctmThreshold = 10; // 10 seconds
    QAction *m_showLogs;
};

class TaskProgressViewPrivate {
public:
    TaskProgressView* q;
    TaskProgressViewModeProxy* proxyModel;
    QToolBar *toolbar;

    QString shortened(QString txt, int maxlen);
};

TaskProgressView::TaskProgressView()
{
    d = new TaskProgressViewPrivate;
    d->q = this;
    d->toolbar = new QToolBar(this);
    setSelectionMode(ContiguousSelection);
    setItemDelegate(new TaskProgressViewDelegate(this));
    TaskManager::TaskProgressModel* model = TaskManager::TaskProgressMonitor::globalInstance()->model();
    d->proxyModel = new TaskProgressViewModeProxy(this);
    d->proxyModel->setSourceModel(model);
    d->proxyModel->setCompleteTasksMode(TaskProgressViewModeProxy::CTM_Hidden);
    setModel(d->proxyModel);
    header()->setSectionResizeMode(0, QHeaderView::Stretch);
    header()->hide();
    setExpandsOnDoubleClick(false);
    connect(this, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(showLogMessages(QModelIndex)));
    QShortcut* copyToClipboard = new QShortcut(QKeySequence(QKeySequence::Copy), this);
    connect(copyToClipboard, SIGNAL(activated()), this, SLOT(copySelectedToClipboard()));
    connect(d->proxyModel, &QAbstractItemModel::modelReset, [this]() {
        for(int i = 0; i < this->model()->rowCount(); ++i)
            this->setFirstColumnSpanned(i, QModelIndex(), true);
    });

    connect(d->proxyModel, &QAbstractItemModel::rowsInserted, [this](const QModelIndex& parent, int from, int to) {
        if (parent.isValid()) return;
        for(int i = from; i <=to; ++i)
            this->setFirstColumnSpanned(i, QModelIndex(), true);
    });
    for (int i = 0; i < this->model()->rowCount(); ++i)
        setFirstColumnSpanned(i, QModelIndex(), true);
    connect(this->model(), &QAbstractItemModel::dataChanged, [this](const QModelIndex& from, const QModelIndex& to) {
        if (from.parent().isValid()) return;
        for (int i = from.row(); i<=to.row(); ++i) {
            setFirstColumnSpanned(i, QModelIndex(), true);
        }
    });
    QTimer* timer = new QTimer(this);
    timer->start(1000);
    connect(timer, SIGNAL(timeout()), viewport(), SLOT(update()));

    // setup actions menu
    setContextMenuPolicy(Qt::ActionsContextMenu);
    QMenu* modeMenu = new QMenu(this);
    modeMenu->setTitle("Mode");
    addAction(modeMenu->menuAction());
    QAction* defMode = new QAction("All tasks", this);
    defMode->setCheckable(true);
    modeMenu->addAction(defMode);
    QAction* activeMode = new QAction("Active tasks", this);
    activeMode->setCheckable(true);
    modeMenu->addAction(activeMode);
//    QAction *lastMinute = new QAction("Active and finished within one minute", this);
//    lastMinute->setCheckable(true);
//    modeMenu->addAction(lastMinute);
    QActionGroup* grp = new QActionGroup(this);
    grp->addAction(defMode);
    grp->addAction(activeMode);
//    grp->addAction(lastMinute);

    d->toolbar->addAction(modeMenu->menuAction());
    QToolButton *modeButton = qobject_cast<QToolButton*>(d->toolbar->widgetForAction(modeMenu->menuAction()));
    if (modeButton) {
        modeButton->setIcon(QIcon::fromTheme("document-preview"));
//    modeButton->setToolTip(modeMenu->title());
        modeButton->setPopupMode(QToolButton::InstantPopup);
    }
    connect(defMode, &QAction::toggled, [this](bool checked) {
       if (checked) {
           d->proxyModel->setCompleteTasksMode(TaskProgressViewModeProxy::CTM_Shown);
       }
    });
    connect(activeMode, &QAction::toggled, [this](bool checked) {
       if (checked) {
           d->proxyModel->setCompleteTasksMode(TaskProgressViewModeProxy::CTM_Hidden);
       }
    });
//    connect(lastMinute, &QAction::toggled, [this](bool checked) {
//       if (checked) {
//           d->proxyModel->setCompleteTasksModeThreshold(60);
//           d->proxyModel->setCompleteTasksMode(TaskProgressViewModeProxy::CTM_HiddenIfOlderThan);
//       }
//    });
    defMode->setChecked(true);
    d->toolbar->addAction(d->proxyModel->showLogsAction());
}

TaskProgressView::~TaskProgressView()
{
    delete d;
}

void TaskProgressView::copySelectedToClipboard()
{
    QItemSelectionModel* selectionModel = this->selectionModel();
    const auto selRows = selectionModel->selectedRows();
    if (selRows.isEmpty())
        return;
    // if first selected entry is a task, we ignore all non-tasks
    bool selectingTasks = !selRows.first().parent().isValid();
    QStringList result;
    QModelIndex lastTask;
    foreach (QModelIndex row, selRows) {
        if (selectingTasks == row.parent().isValid())
            continue;
        if (selectingTasks) {
            // for each task get the name of the task
            result << row.data().toString();
            // for each task get its log messages
            result << row.data(TaskManager::TaskProgressModel::IndentedLogRole).toString();
        }
        else {
            // for each log see if it belongs to the previos task
            // if not, add the task name to the log
            if (row.parent() != lastTask) {
                result << row.parent().data().toString();
                lastTask = row.parent();
            }
            result << row.data(TaskManager::TaskProgressModel::IndentedLogRole).toString();
        }
    }
    QApplication::clipboard()->setText(result.join("\n"));
}

void TaskProgressView::showLogMessages(const QModelIndex& index)
{
    if (index.parent().isValid())
        return;
    QWidget* dlg = new QWidget(this);
    // forcing the widget to be a window despite having a parent so that the window
    // is destroyed with the main window
    dlg->setWindowFlags(Qt::Window);
    dlg->setAttribute(Qt::WA_DeleteOnClose);
    dlg->setWindowTitle(tr("Log messages for %1").arg(index.data().toString()));
    QPlainTextEdit* te = new QPlainTextEdit;
    QFont f;
    f.setPointSize(f.pointSize() - 2);
    te->setFont(f);
    te->setReadOnly(true);
    QDialogButtonBox* bb = new QDialogButtonBox;
    QObject::connect(bb, SIGNAL(accepted()), dlg, SLOT(close()));
    QObject::connect(bb, SIGNAL(rejected()), dlg, SLOT(close()));
    bb->setStandardButtons(QDialogButtonBox::Close);
    QVBoxLayout* l = new QVBoxLayout(dlg);
    l->addWidget(te);
    l->addWidget(bb);
    te->setPlainText(index.data(TaskManager::TaskProgressModel::LogRole).toString());
    QRect r = QApplication::desktop()->screenGeometry(dlg);
    r.setWidth(r.width() / 2);
    r.setHeight(r.height() / 2);
    dlg->resize(r.size());
    dlg->show();
}

void TaskProgressView::updateGeometries()
{
    QSize tbSize = d->toolbar->sizeHint();
    d->toolbar->setGeometry(0, 0, width(), tbSize.height());
    setViewportMargins(0, tbSize.height(), 0, 0);
    QTreeView::updateGeometries();
    QRect scrollbarRect = verticalScrollBar()->geometry();
    scrollbarRect.setTop(tbSize.height());
    verticalScrollBar()->setGeometry(scrollbarRect);
}

QString TaskProgressViewPrivate::shortened(QString txt, int maxlen)
{
    if (txt.count() > maxlen) {
        return txt.mid(0, maxlen - 3) + "...";
    }
    else
        return txt;
}
