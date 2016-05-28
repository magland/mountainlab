/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/30/2016
*******************************************************/

#include "taskprogress.h"

#include <QTimer>
#include <QCoreApplication>
#include <QDebug>
#include <QMutex>

class TaskProgressPrivate {
public:
    TaskProgress* q;
    TaskInfo m_info;
};

class TaskProgressAgentPrivate {
public:
    TaskProgressAgent* q;
    QList<TaskInfo> m_completed_tasks;
    QList<TaskProgress*> m_active_tasks;
    bool m_emit_tasks_changed_scheduled;
    QDateTime m_last_emit_tasks_changed;
    QMutex m_addtask_removetask_mutex;
    QMap<QString, double> m_quantities;
    QMutex m_mutex;
};

TaskProgress::TaskProgress(const QString& label, const QString& description)
{
    d = new TaskProgressPrivate;
    d->q = this;

    d->m_info.label = label;
    d->m_info.description = description;
    d->m_info.progress = 0;
    d->m_info.start_time = QDateTime::currentDateTime();

    TaskProgressAgent::globalInstance()->addTask(this);
}

TaskProgress::~TaskProgress()
{
    this->log("Destructor");
    if (d->m_info.progress != 1) {
        d->m_info.end_time = QDateTime::currentDateTime();
    }
    TaskProgressAgent::globalInstance()->removeTask(this);
    delete d;
}

void TaskProgress::setLabel(const QString& label)
{
    {
        if (d->m_info.label == label)
            return;
        d->m_info.label = label;
        emit changed();
    }
    this->log(label);
}

void TaskProgress::setDescription(const QString& description)
{
    {
        if (d->m_info.description == description)
            return;
        d->m_info.description = description;
        emit changed();
    }
    this->log(description);
}

void TaskProgress::log(const QString& log_message)
{
    TaskProgressLogMessage MSG;
    MSG.message = log_message;
    MSG.time = QDateTime::currentDateTime();
    d->m_info.log_messages << MSG;
    emit changed();
}

void TaskProgress::error(const QString& error_message)
{
    qWarning() << "TaskProgress Error:: " + error_message;
    this->log("ERROR: " + error_message);
    {
        d->m_info.error = error_message;
    }
}

void TaskProgress::setProgress(double pct)
{
    if (d->m_info.progress == pct)
        return;
    d->m_info.progress = pct;
    emit changed();
    if (pct == 1) {
        d->m_info.end_time = QDateTime::currentDateTime();
        emit completed(d->m_info);
    }
}

TaskInfo TaskProgress::getInfo() const
{
    return d->m_info;
}

TaskProgressAgent::TaskProgressAgent()
{
    d = new TaskProgressAgentPrivate;
    d->q = this;
    d->m_emit_tasks_changed_scheduled = false;
}

TaskProgressAgent::~TaskProgressAgent()
{
    delete d;
}

QList<TaskInfo> TaskProgressAgent::activeTasks()
{
    QMutexLocker locker(&d->m_mutex);
    QList<TaskInfo> ret;
    QList<TaskProgress*> to_remove;
    foreach(TaskProgress * X, d->m_active_tasks)
    {
        TaskInfo info = X->getInfo();
        if (!info.label.isEmpty()) {
            if (info.progress < 1) {
                ret << X->getInfo();
            }
        }
    }

    return ret;
}

QList<TaskInfo> TaskProgressAgent::completedTasks()
{
    QMutexLocker locker(&d->m_mutex);
    for (int i = 0; i < d->m_completed_tasks.count(); i++) {
        d->m_completed_tasks[i].progress = 1; // kind of a hack to make sure the progress is 1 for all completed tasks
    }
    QList<TaskInfo> ret;
    for (int i = 0; i < d->m_completed_tasks.count(); i++) {
        if (!d->m_completed_tasks[i].label.isEmpty())
            ret << d->m_completed_tasks[i];
    }
    return ret;
}

/// Witold will this way work?
Q_GLOBAL_STATIC(QMutex, global_instance_mutex)
Q_GLOBAL_STATIC(TaskProgressAgent, theInstance)
TaskProgressAgent* TaskProgressAgent::globalInstance()
{
    QMutexLocker locker(global_instance_mutex);
    return theInstance;
}

void TaskProgressAgent::slot_schedule_emit_tasks_changed()
{
    if (d->m_emit_tasks_changed_scheduled)
        return;
    if (d->m_last_emit_tasks_changed.secsTo(QDateTime::currentDateTime()) >= 1) {
        slot_emit_tasks_changed();
    } else {
        d->m_emit_tasks_changed_scheduled = true;
        QTimer::singleShot(300, this, SLOT(slot_emit_tasks_changed()));
    }
}

void TaskProgressAgent::slot_emit_tasks_changed()
{
    {
        d->m_last_emit_tasks_changed = QDateTime::currentDateTime();
        d->m_emit_tasks_changed_scheduled = false;
    }
    emit tasksChanged();
}

void TaskProgressAgent::slot_task_completed(TaskInfo info)
{
    QMutexLocker locker(&d->m_mutex);
    d->m_completed_tasks.prepend(info);
}

void TaskProgressAgent::addTask(TaskProgress* T)
{
    {
        QMutexLocker locker(&d->m_mutex);
        QString label = T->getInfo().label;
        d->m_active_tasks.prepend(T);
        connect(T, SIGNAL(changed()), this, SLOT(slot_schedule_emit_tasks_changed()), Qt::QueuedConnection);
        connect(T, SIGNAL(completed(TaskInfo)), this, SLOT(slot_task_completed(TaskInfo)), Qt::QueuedConnection);
    }
    emit tasksChanged();
}

void TaskProgressAgent::removeTask(TaskProgress* T)
{
    {
        QMutexLocker locker(&d->m_mutex);
        QString label = T->getInfo().label;
        TaskInfo info = T->getInfo();
        if (info.progress != 1) {
            d->m_completed_tasks.prepend(T->getInfo());
        }
        d->m_active_tasks.removeAll(T);
    }
    emit tasksChanged();
}

void TaskProgressAgent::incrementQuantity(QString name, double val)
{
    {
        QMutexLocker locker(&d->m_mutex);
        d->m_quantities[name] = d->m_quantities[name] + val;
    }
    emit quantitiesChanged();
}

double TaskProgressAgent::getQuantity(QString name)
{
    QMutexLocker locker(&d->m_mutex);
    return d->m_quantities.value(name);
}
