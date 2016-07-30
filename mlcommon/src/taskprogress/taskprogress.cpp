/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland, Witold Wysota
** Created: 4/30/2016
*******************************************************/

#include "taskprogress.h"

#include <QTimer>
#include <QCoreApplication>
#include <QDebug>
#include <QMutex>
#include <QThread>
#include <QSortFilterProxyModel>

#include <QAtomicInt>

Q_DECLARE_METATYPE(QSet<QString>)
namespace TaskManager {

class TaskProgressAgentPrivate : public TaskProgressAgent {
public:
    TaskProgressAgentPrivate(int id, const TaskInfo& info)
        : TaskProgressAgent()
        , m_id(id)
        , m_info(info)
    {
    }
    const TaskInfo& taskInfo() const override
    {
        return m_info;
    }
    void emitChanged() { emit changed(); }

    int m_id;
    TaskInfo m_info;
};

class TaskProgressModelPrivate : public TaskProgressModel {
public:
    using TaskProgressModel::index;

    TaskProgressModelPrivate(QObject* parent = 0)
        : TaskProgressModel(parent)
    {
    }

    void addTask(int id, TaskInfo info)
    {
        int from = 0;
        beginInsertRows(QModelIndex(), from, from);
        info.start_time = QDateTime::currentDateTime();
        TaskProgressAgentPrivate* nInfo = new TaskProgressAgentPrivate(id, info);
        m_map.insert(id, nInfo);
        m_data.insert(from, nInfo);
        endInsertRows();
        m_activeCount++;
    }

    void modifyTask(const QModelIndex& task, const QVariant& value, int role = Qt::EditRole)
    {
        if (!isTask(task))
            return;
        int row = task.row();
        TaskProgressAgentPrivate* info = m_data.at(task.row());
        switch (role) {
        case TaskProgressModel::ProgressRole:
            info->m_info.progress = value.toDouble();
            if (info->m_info.progress >= 1) {
                completeTask(task);
            } else {
                info->emitChanged();
        emit dataChanged(index(row, 0), index(row, columnCount() - 1));
            }
        }
    }

    void setLabel(const QModelIndex& idx, const QString& label)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        if (task->m_info.label == label)
            return;
        task->m_info.label = label;
        task->emitChanged();
        emit dataChanged(index(idx.row(), 0), index(idx.row(), columnCount() - 1));
    }

    void setDescription(const QModelIndex& idx, const QString& description)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        if (task->m_info.description == description)
            return;
        task->m_info.description = description;
        task->emitChanged();
        emit dataChanged(index(idx.row(), 0), index(idx.row(), columnCount() - 1));
    }

    void log(const QModelIndex& idx, const QString& message)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        beginInsertRows(idx, 0, 0);
        task->m_info.log_messages.append(message);
        endInsertRows();
        task->emitChanged();
    }

    void error(const QModelIndex& idx, const QString& message)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        beginInsertRows(idx, 0, 0);
        task->m_info.log_messages.append(tr("Error: %1").arg(message));
        task->m_info.error = message;
        endInsertRows();
        task->emitChanged();
    }

    void addTag(const QModelIndex& idx, const QString& tag)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        task->m_info.tags << tag;
        task->emitChanged();
        emit dataChanged(index(idx.row(), 0), index(idx.row(), columnCount() - 1));
    }

    void removeTag(const QModelIndex& idx, const QString& tag)
    {
        TaskProgressAgentPrivate* task = info(idx);
        if (!task)
            return;
        task->m_info.tags.remove(tag);
        task->emitChanged();
        emit dataChanged(index(idx.row(), 0), index(idx.row(), columnCount() - 1));
    }

    void completeTask(const QModelIndex& index)
    {
        TaskProgressAgentPrivate* task = info(index);
        if (!task)
            return;
        int row = index.row();
        int newRow = rowCount();
        task->m_info.progress = 1.0;
        task->m_info.end_time = QDateTime::currentDateTime();
        task->emitChanged();
        if (row == newRow - 1)
        {
            emit dataChanged(this->index(row, 0), this->index(row, columnCount() - 1));
            return;
        }
        beginMoveRows(QModelIndex(), row, row, QModelIndex(), newRow);
        m_data.move(row, newRow - 1);
        m_activeCount--;
        endMoveRows();
    }

    TaskProgressAgentPrivate* info(const QModelIndex& index) const
    {
        if (isTask(index))
            return m_data.at(index.row());
        return 0;
    }
    QModelIndex index(TaskProgressAgent* agent) const
    {
        int idx = m_data.indexOf(static_cast<TaskProgressAgentPrivate*>(agent));
        if (idx < 0)
            return QModelIndex();
        return index(idx, 0);
    }

    TaskProgressAgentPrivate* infoById(int id) const
    {
        return m_map.value(id, 0);
    }

    QModelIndex indexById(int id) const
    {
        TaskProgressAgentPrivate* agent = infoById(id);
        return index(agent);
    }

    void cancelTask(const QModelIndex& index)
    {
        completeTask(index);
    }

    QHash<int, TaskProgressAgentPrivate*> m_map;
};

class TaskProgressEvent : public QEvent {
public:
    enum Command {
        Create,
        SetLabel,
        SetDescription,
        AddTag,
        RemoveTag,
        AppendLog,
        AppendError,
        SetProgress,
        Finish
    };
    static QEvent::Type type()
    {
        static QEvent::Type typeVal = static_cast<QEvent::Type>(registerEventType());
        return typeVal;
    }

    TaskProgressEvent(int id, Command cmd, const QVariant& value = QVariant())
        : QEvent(TaskProgressEvent::type())
    {
        m_id = id;
        m_cmd = cmd;
        m_value = value;
    }
    Command command() const { return m_cmd; }
    const QVariant& value() const { return m_value; }
    int id() const { return m_id; }

private:
    int m_id;
    Command m_cmd;
    QVariant m_value;
};

class TaskProgressMonitorPrivate : public TaskProgressMonitor {
public:
    TaskProgressMonitorPrivate()
    {
        m_model = new TaskProgressModelPrivate(this);
        qRegisterMetaType<TaskInfo>();
    }

    ~TaskProgressMonitorPrivate()
    {
        delete m_model;
    }

    void post(TaskProgressEvent* e)
    {
        QCoreApplication::postEvent(this, e);
    }

    int count() const override
    {
        return m_model->rowCount();
    }

    TaskProgressAgent* at(int index) const override
    {
        return m_model->info(m_model->index(index, 0));
    }

    int indexOf(TaskProgressAgent* agent) const override
    {
        QModelIndex idx = m_model->index(agent);
        if (idx.isValid())
            return idx.row();
        return -1;
    }

    void incrementQuantity(QString name, double val) override
    {
        m_quantities.insert(name, m_quantities.value(name, 0) + val);
        emit quantitiesChanged();
    }
    double getQuantity(QString name) const override
    {
        return m_quantities.value(name, 0);
    }

    TaskProgressModel* model() const override
    {
        return m_model;
    }

    void complete(int idx)
    {
        m_model->completeTask(m_model->index(idx, 0));
    }

    static TaskProgressMonitorPrivate* privateInstance();

protected:
    void customEvent(QEvent* event)
    {
        if (event->type() != TaskProgressEvent::type()) {
            event->ignore();
            return;
        }
        TaskProgressEvent* tpEvent = static_cast<TaskProgressEvent*>(event);
        if (tpEvent->command() == TaskProgressEvent::Create) {
            TaskInfo info = tpEvent->value().value<TaskInfo>();
            m_model->addTask(tpEvent->id(), info);
        }
        QModelIndex index = m_model->indexById(tpEvent->id());
        const QVariant& value = tpEvent->value();
        switch (tpEvent->command()) {
        case TaskProgressEvent::SetLabel:
            m_model->setLabel(index, value.toString());
            break;
        case TaskProgressEvent::SetDescription:
            m_model->setDescription(index, value.toString());
            break;
        case TaskProgressEvent::AddTag:
            m_model->addTag(index, value.toString());
            break;
        case TaskProgressEvent::RemoveTag:
            m_model->removeTag(index, value.toString());
            break;
        case TaskProgressEvent::AppendLog:
            m_model->log(index, value.toString());
            break;
        case TaskProgressEvent::AppendError:
            m_model->error(index, value.toString());
            break;
        case TaskProgressEvent::SetProgress:
            m_model->modifyTask(index, value, TaskProgressModel::ProgressRole);
            break;
        case TaskProgressEvent::Finish:
            m_model->completeTask(index);
            break;
        case TaskProgressEvent::Create:
            break; // handled earlier
        default:
            qWarning("Custom event not handled");
            break;
        }
    }

private:
    TaskProgressModelPrivate* m_model;
    QMap<QString, double> m_quantities;
};

Q_GLOBAL_STATIC(TaskProgressMonitorPrivate, _q_tpm_instance)

TaskProgressMonitor* TaskProgressMonitor::globalInstance()
{
    return _q_tpm_instance;
}

TaskProgressMonitorPrivate* TaskProgressMonitorPrivate::privateInstance()
{
    return _q_tpm_instance;
}

class TagsFilterProxyModel : public QSortFilterProxyModel {
public:
    TagsFilterProxyModel(QObject* parent = 0)
        : QSortFilterProxyModel(parent)
    {
    }
    void addTag(const QString& t)
    {
        m_tags << t;
        invalidateFilter();
    }
    void setTags(const QSet<QString>& tags)
    {
        m_tags = tags;
        invalidateFilter();
    }

protected:
    bool filterAcceptsRow(int source_row, const QModelIndex& source_parent) const
    {
        if (m_tags.isEmpty())
            return true;
        if (source_parent.isValid())
            return true;
        TaskProgressAgent* a = TaskProgressMonitor::globalInstance()->at(source_row);
#if QT_VERSION >= 0x050600
        return m_tags.intersects(a->taskInfo().tags);
#else
        return !QSet<QString>(m_tags).intersect(a->taskInfo().tags).isEmpty();
#endif
    }

private:
    QSet<QString> m_tags;
};

TaskProgressMonitorFilter::TaskProgressMonitorFilter(QObject* parent)
    : QObject(parent)
{
    m_proxy = new TagsFilterProxyModel(this);
    m_base = TaskProgressMonitor::globalInstance()->model();
    m_proxy->setSourceModel(m_base);
}

TaskProgressMonitorFilter::TaskProgressMonitorFilter(const QString& tag, QObject* parent)
    : QObject(parent)
{
    m_proxy = new TagsFilterProxyModel(this);
    m_proxy->addTag(tag);
    m_base = TaskProgressMonitor::globalInstance()->model();
    m_proxy->setSourceModel(m_base);
}

TaskProgressMonitorFilter::TaskProgressMonitorFilter(const QStringList& tags, QObject* parent)
    : QObject(parent)
{
    m_proxy = new TagsFilterProxyModel(this);
    m_proxy->setTags(tags.toSet());
    m_base = TaskProgressMonitor::globalInstance()->model();
    m_proxy->setSourceModel(m_base);
}

void TaskProgressMonitorFilter::addTag(const QString& tag)
{
    m_proxy->addTag(tag);
}

int TaskProgressMonitorFilter::count() const
{
    return model()->rowCount();
}

TaskProgressAgent* TaskProgressMonitorFilter::at(int index) const
{
    QModelIndex idx = model()->index(index, 0);
    QModelIndex srcIdx = m_proxy->mapToSource(idx);
    if (!srcIdx.isValid())
        return 0;
    return TaskProgressMonitor::globalInstance()->at(srcIdx.row());
}

int TaskProgressMonitorFilter::indexOf(TaskProgressAgent* a) const
{
    int srcRow = TaskProgressMonitor::globalInstance()->indexOf(a);
    if (srcRow < 0)
        return -1;
    QModelIndex srcIdx = m_proxy->mapFromSource(m_base->index(srcRow, 0));
    if (srcIdx.isValid())
        return srcIdx.row();
    return -1;
}

QAbstractItemModel* TaskProgressMonitorFilter::model() const
{
    return m_proxy;
}

TaskProgressModel::TaskProgressModel(QObject* parent)
    : QAbstractItemModel(parent)
{
}

TaskProgressModel::~TaskProgressModel()
{
    qDeleteAll(m_data);
}

QModelIndex TaskProgressModel::index(int row, int column, const QModelIndex& parent) const
{
    if (!parent.isValid()) {
        // task
        if (row < 0 || row >= m_data.count())
            return QModelIndex();
        return createIndex(row, column, nullptr);
    }
    if (parent.internalPointer() == 0) {
        TaskProgressAgentPrivate* info = m_data.at(parent.row());
        if (row >= info->taskInfo().log_messages.count() || row < 0)
            return QModelIndex();
        // log entry
        return createIndex(row, column, info);
    }
    qFatal("Should never reach this line");
    return QModelIndex();
}

QModelIndex TaskProgressModel::parent(const QModelIndex& child) const
{
    if (!child.isValid())
        return QModelIndex();
    if (child.internalPointer() == 0)
        return QModelIndex();
    TaskProgressAgentPrivate* info = static_cast<TaskProgressAgentPrivate*>(child.internalPointer());
    return createIndex(m_data.indexOf(info), 0, nullptr);
}

int TaskProgressModel::rowCount(const QModelIndex& parent) const
{
    if (!parent.isValid())
        return m_data.size();
    if (!parent.internalPointer()) {
        TaskProgressAgentPrivate* info = m_data.at(parent.row());
        return info->taskInfo().log_messages.size();
    }
    return 0;
}

int TaskProgressModel::columnCount(const QModelIndex&) const
{
    return 2;
}

QVariant TaskProgressModel::data(const QModelIndex& index, int role) const
{
    if (!index.isValid())
        return QVariant();
    if (isTask(index))
        return taskData(index, role);
    return logData(index, role);
}

QVariant TaskProgressModel::taskData(const QModelIndex& index, int role) const
{
    if (index.column() != 0)
        return QVariant();
    const TaskInfo& task = m_data.at(index.row())->taskInfo();
    switch (role) {
    case Qt::EditRole:
    case Qt::DisplayRole:
        // modified by jfm -- 5/17/2016
        if (!task.error.isEmpty()) {
            return QString("%1: %2").arg(task.label).arg(task.error);
        }
        else {
            return task.label;
        }
    case Qt::ToolTipRole:
        // modified by jfm -- 5/17/2016
        if (!task.description.isEmpty())
            return task.description;
        else
            return taskData(index, Qt::DisplayRole);
    ////#ifdef QT_WIDGETS_LIB
    //    case Qt::ForegroundRole: {
    //        // modified by jfm -- 5/17/2016
    //        if (!task.error.isEmpty()) {
    //            return QColor(Qt::red);
    //        }
    //        else {
    //            if (task.progress < 1)
    //                return QColor(Qt::blue);
    //        }
    //        return QVariant();
    //    }
    ////#endif // mlcommon doesn't link QtGui anyway
    case ProgressRole:
        return task.progress;
    case StartTimeRole:
        return task.start_time;
    case EndTimeRole:
        return task.end_time;
    case TagsRole:
        return QStringList(task.tags.toList());
    case TagsSetRole:
        return qVariantFromValue(task.tags);
    case LogRole:
        return assembleLog(task);
    case IndentedLogRole:
        return assembleLog(task, "\t");
    case StatusRole:
        if (task.end_time.isValid() && task.progress >= 1)
            return Finished;
        if (task.end_time.isValid())
            return Canceled;
        return Active;
    }
    return QVariant();
}

QVariant TaskProgressModel::logData(const QModelIndex& index, int role) const
{

    TaskProgressAgentPrivate* agent = static_cast<TaskProgressAgentPrivate*>(index.internalPointer());
    const auto& logMessages = agent->taskInfo().log_messages;
    int idx = logMessages.count() - 1 - index.row();

    auto logMessage = logMessages.at(idx); // newest first
    switch (role) {
    case Qt::EditRole:
    case Qt::DisplayRole:
        //        return logMessage.time.toString(Qt::SystemLocaleShortDate) + "\t" + logMessage.message;
        if (index.column() == 0)
            return logMessage.time;
        return logMessage.message;
    case Qt::UserRole:
        return logMessage.time;
    case LogRole:
        return singleLog(logMessage);
    case IndentedLogRole:
        return singleLog(logMessage, "\t");
    default:
        return QVariant();
    }
}

bool TaskProgressModel::isActive(const QModelIndex& task) const
{
    if (!isTask(task))
        return false;
    return (task.data(StatusRole).toInt() == Active);
}

bool TaskProgressModel::isCompletedWithin(const QModelIndex& task, int time) const
{
    if (!isTask(task))
        return false;
    Status s = (Status)task.data(StatusRole).toInt();
    if (s == Active)
        return true;
    const QDateTime dt = task.data(EndTimeRole).toDateTime();
    if (dt.addSecs(time) >= QDateTime::currentDateTime())
        return true;
    return false;
}

bool TaskProgressModel::isTask(const QModelIndex& idx) const
{
    return (idx.isValid() && idx.internalPointer() == nullptr);
}

QString TaskProgressModel::assembleLog(const TaskInfo& task, const QString& prefix) const
{
    QStringList entries;
    foreach (const TaskProgressLogMessage& msg, task.log_messages) {
        entries << singleLog(msg, prefix);
    }
    return entries.join("\n");
}

QString TaskProgressModel::singleLog(const TaskProgressLogMessage& msg, const QString& prefix) const
{
    return QString("%1%2: %3").arg(prefix).arg(msg.time.toString("yyyy-MM-dd hh:mm:ss.zzz")).arg(msg.message);
}
}

static QAtomicInt TaskProgressValue = 0;

TaskProgress::TaskProgress()
    : QObject()
    , m_id(TaskProgressValue.fetchAndAddOrdered(1))
{
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::Create);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

TaskProgress::TaskProgress(const QString& label)
    : QObject()
    , m_id(TaskProgressValue.fetchAndAddOrdered(1))
{
    m_info.label = label;
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::Create, qVariantFromValue(m_info));
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

TaskProgress::TaskProgress(StandardCategories tags, const QString& label)
    : QObject()
    , m_id(TaskProgressValue.fetchAndAddOrdered(1))
{
    m_info.label = label;
    QStringList tagNames = catsToString(tags);
    m_info.tags = tagNames.toSet();
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::Create, qVariantFromValue(m_info));
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

TaskProgress::~TaskProgress()
{
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::Finish);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

QString TaskProgress::label() const
{
    return m_info.label;
}

QString TaskProgress::description() const
{
    return m_info.description;
}

QSet<QString> TaskProgress::tags() const
{
    return m_info.tags;
}

void TaskProgress::setLabel(const QString& label)
{
    if (m_info.label == label)
        return;
    m_info.label = label;
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::SetLabel,
        label);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

void TaskProgress::setDescription(const QString& description)
{
    if (m_info.description == description)
        return;
    m_info.description = description;
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::SetDescription,
        description);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

void TaskProgress::addTag(TaskProgress::StandardCategory cat)
{
    QString catName = catToString(cat);
    if (catName.isEmpty())
        return;
    m_info.tags << catName;
}

void TaskProgress::addTag(const QString& tag)
{
    m_info.tags << tag;
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::AddTag,
        tag);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

void TaskProgress::removeTag(TaskProgress::StandardCategory cat)
{
    QString catName = catToString(cat);
    if (catName.isEmpty())
        return;
    removeTag(catName);
}

void TaskProgress::removeTag(const QString& tag)
{
    m_info.tags.remove(tag);
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::RemoveTag,
        tag);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

bool TaskProgress::hasTag(TaskProgress::StandardCategory cat) const
{
    QString catName = catToString(cat);
    if (catName.isEmpty())
        return false;
    return hasTag(catName);
}

bool TaskProgress::hasTag(const QString& tag) const
{
    return m_info.tags.contains(tag);
}

double TaskProgress::progress() const
{
    return m_info.progress;
}

void TaskProgress::log(const QString& log_message)
{
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::AppendLog,
        log_message);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

void TaskProgress::error(const QString& error_message)
{
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::AppendError,
        error_message);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

void TaskProgress::setProgress(double pct)
{
    if (m_info.progress == pct)
        return;
    m_info.progress = pct;
    TaskManager::TaskProgressEvent* event = new TaskManager::TaskProgressEvent(m_id,
        TaskManager::TaskProgressEvent::SetProgress,
        pct);
    TaskManager::TaskProgressMonitorPrivate::privateInstance()->post(event);
}

QStringList TaskProgress::catsToString(StandardCategories cats) const
{
    QStringList result;
    if (cats.testFlag(Download))
        result << "download";
    if (cats.testFlag(Calculate))
        result << "calculate";
    if (cats.testFlag(Process))
        result << "process";
    return result;
}

QString TaskProgress::catToString(TaskProgress::StandardCategory cat) const
{
    switch (cat) {
    case Download:
        return "download";
    case Calculate:
        return "calculate";
    case Process:
        return "process";
    case None:
    default:
        return QString();
    }
}
