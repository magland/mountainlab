/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/30/2016
*******************************************************/

#ifndef TASKPROGRESS_H
#define TASKPROGRESS_H

#include <QDateTime>
#include <QObject>
#include <QString>
#include <QDebug>

struct TaskProgressLogMessage {
    QString message;
    QDateTime time;
};

struct TaskInfo {
    QString label;
    QString description;
    QList<TaskProgressLogMessage> log_messages;
    QString error;
    double progress;
    QDateTime start_time;
    QDateTime end_time;
};

Q_DECLARE_METATYPE(TaskInfo)

class TaskProgressPrivate;
class TaskProgress : public QObject {
    Q_OBJECT
public:
    friend class TaskProgressPrivate;
    TaskProgress(const QString& label = "", const QString& description = "");
    virtual ~TaskProgress();
    void setLabel(const QString& label);
    void setDescription(const QString& description);
    void log(const QString& log_message);
    void error(const QString& error_message);
    void setProgress(double pct);

    TaskInfo getInfo() const;
signals:
    void changed();
    void completed(TaskInfo info);

private:
    TaskProgressPrivate* d;
};

class TaskProgressAgentPrivate;
class TaskProgressAgent : public QObject {
    Q_OBJECT
public:
    friend class TaskProgressAgentPrivate;
    TaskProgressAgent();
    virtual ~TaskProgressAgent();
    void addTask(TaskProgress* T);
    void removeTask(TaskProgress* T);
    void incrementQuantity(QString name, double val);
    double getQuantity(QString name);

    QList<TaskInfo> activeTasks();
    QList<TaskInfo> completedTasks();

    static TaskProgressAgent* globalInstance();
signals:
    void tasksChanged();
    void quantitiesChanged();
private slots:
    void slot_schedule_emit_tasks_changed();
    void slot_emit_tasks_changed();
    void slot_task_completed(TaskInfo info);

private:
    TaskProgressAgentPrivate* d;
};

#endif // TASKPROGRESS_H
