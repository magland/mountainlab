#include "prvmanagerdialog.h"
#include "ui_prvmanagerdialog.h"

#include <QProcess>
#include <QDebug>
#include <QThread>
#include <QMutex>
#include <taskprogress.h>
#include <QMessageBox>
#include <QInputDialog>
#include <cachemanager.h>
#include <QJsonDocument>
#include "mlcommon.h"

class PrvManagerDialogPrivate {
public:
    PrvManagerDialog* q;
    QMap<QString, QJsonObject> m_prv_objects;
    QStringList m_server_names;
    PrvManagerDialogThread m_thread;

    void refresh_tree();
    void ensure_local(QJsonObject prv_obj);
    void ensure_remote(QJsonObject prv_obj, QString server);

    void restart_thread();
};

PrvManagerDialog::PrvManagerDialog(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::PrvManagerDialog)
{
    d = new PrvManagerDialogPrivate;
    d->q = this;
    ui->setupUi(this);

    QObject::connect(ui->btn_upload_to_server, SIGNAL(clicked(bool)), this, SLOT(slot_upload_to_server()));
    QObject::connect(ui->btn_copy_to_local_disk, SIGNAL(clicked(bool)), this, SLOT(slot_copy_to_local_disk()));
    QObject::connect(ui->btn_done, SIGNAL(clicked(bool)), this, SLOT(close()));
    QObject::connect(&d->m_thread, SIGNAL(results_updated()), this, SLOT(slot_results_updated()));
}

PrvManagerDialog::~PrvManagerDialog()
{
    if (d->m_thread.isRunning()) {
        d->m_thread.requestInterruption();
        d->m_thread.wait(5000);
        d->m_thread.terminate();
    }
    delete ui;
    delete d;
}

void PrvManagerDialog::setPrvObjects(const QMap<QString, QJsonObject>& prv_objects)
{
    d->m_prv_objects = prv_objects;
}

void PrvManagerDialog::setServerNames(const QStringList& server_names)
{
    d->m_server_names = server_names;
}

void PrvManagerDialog::refresh()
{
    d->refresh_tree();
}

QString to_string(fuzzybool fb)
{
    if (fb == YES)
        return "YES";
    if (fb == NO)
        return "x";
    return ".";
}

void PrvManagerDialog::slot_results_updated()
{
    d->m_thread.results_mutex.lock();
    QMap<QString, PrvManagerDialogResult> results = d->m_thread.results;
    d->m_thread.results_mutex.unlock();

    QTreeWidget* TT = ui->tree;
    for (int i = 0; i < TT->topLevelItemCount(); i++) {
        QTreeWidgetItem* it = TT->topLevelItem(i);
        QString name = it->data(0, Qt::UserRole).toString();

        PrvManagerDialogResult result0 = results[name];
        {
            int col = 2;
            it->setText(col, to_string(result0.on_local_disk));
        }
        for (int a = 0; a < d->m_server_names.count(); a++) {
            int col = 3 + a;
            it->setText(col, to_string(result0.on_server.value(d->m_server_names[a], UNKNOWN)));
        }
    }
}

void PrvManagerDialog::slot_upload_to_server()
{
    QTreeWidget* TT = ui->tree;
    QTreeWidgetItem* it = TT->currentItem();
    if (!it) {
        QMessageBox::warning(0, "Upload to server", "You must select an item first.");
        return;
    }
    QString server = QInputDialog::getText(0, "Upload to server", "Name of server:");
    if (server.isEmpty())
        return;
    QString name = it->data(0, Qt::UserRole).toString();
    QJsonObject obj = d->m_prv_objects[name];
    d->ensure_remote(obj, server);
}

void PrvManagerDialog::slot_copy_to_local_disk()
{
    QTreeWidget* TT = ui->tree;
    QTreeWidgetItem* it = TT->currentItem();
    if (!it) {
        QMessageBox::warning(0, "Copy to local disk", "You must select an item first.");
        return;
    }
    QString name = it->data(0, Qt::UserRole).toString();
    QJsonObject obj = d->m_prv_objects[name];
    d->ensure_local(obj);
}

void PrvManagerDialog::slot_restart_thread()
{
    d->restart_thread();
}

QString format_file_size(int file_size)
{
    if (file_size < 1e3) {
        return QString("%1 bytes").arg(file_size);
    }
    else if (file_size < 1e6) {
        return QString("%1K").arg(((int)(file_size * 1.0 / 1e3 * 10)) / 10.0);
    }
    else if (file_size < 1e9) {
        return QString("%1M").arg(((int)(file_size * 1.0 / 1e6 * 10)) / 10.0);
    }
    else {
        return QString("%1G").arg(((int)(file_size * 1.0 / 1e9 * 10)) / 10.0);
    }
}

void PrvManagerDialogPrivate::refresh_tree()
{
    QTreeWidget* TT = q->ui->tree;
    TT->clear();
    QStringList names = m_prv_objects.keys();

    QStringList headers;
    headers << "Name"
            << "Size"
            << "Local disk";
    headers.append(m_server_names);
    TT->setHeaderLabels(headers);

    for (int i = 0; i < names.count(); i++) {
        QString name = names[i];
        QJsonObject prv_object = m_prv_objects[name];
        bigint file_size = prv_object["original_size"].toVariant().toLongLong();
        QTreeWidgetItem* it = new QTreeWidgetItem();
        it->setText(0, name);
        it->setData(0, Qt::UserRole, name);
        it->setText(1, format_file_size(file_size));
        for (int col = 1; col < headers.count(); col++) {
            it->setTextAlignment(col, Qt::AlignCenter);
        }
        TT->addTopLevelItem(it);
    }

    for (int i = 0; i < TT->columnCount(); i++) {
        TT->resizeColumnToContents(i);
    }

    restart_thread();
}

class ExecCmdThread : public QThread {
public:
    QString cmd;
    QStringList args;

    void run()
    {
        TaskProgress task("Running: " + cmd + " " + args.join(" "));
        task.log() << "Running: " + cmd + " " + args.join(" ");
        QProcess P;
        P.setReadChannelMode(QProcess::MergedChannels);
        P.start(cmd, args);
        P.waitForStarted();
        P.waitForFinished(-1);
        task.log() << "Exit code: " + P.exitCode();
        task.log() << P.readAll();
    }
};

void execute_command_in_separate_thread(QString cmd, QStringList args, QObject* on_finished_receiver, const char* sig_or_slot)
{
    ExecCmdThread* thread = new ExecCmdThread;

    if (on_finished_receiver) {
        QObject::connect(thread, SIGNAL(finished()), on_finished_receiver, sig_or_slot);
    }

    /// Witold, is the following line okay?
    QObject::connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));

    thread->cmd = cmd;
    thread->args = args;
    thread->start();
}

void PrvManagerDialogPrivate::ensure_local(QJsonObject prv_obj)
{
    QString tmp_fname = CacheManager::globalInstance()->makeLocalFile(MLUtil::makeRandomId(10) + ".PrvManagerdlg.prv");
    TextFile::write(tmp_fname, QJsonDocument(prv_obj).toJson());
    CacheManager::globalInstance()->setTemporaryFileDuration(tmp_fname, 600);
    QString cmd = "prv";
    QStringList args;
    args << "ensure-local" << tmp_fname;
    execute_command_in_separate_thread(cmd, args, q, SLOT(slot_restart_thread()));

    /*
    if (QProcess::execute(cmd, args) == 0) {
        restart_thread();
    }
    else {
        QMessageBox::warning(0,"Ensure local","Problem running: "+cmd+" "+args.join(" "));
    }
    */
}

void PrvManagerDialogPrivate::ensure_remote(QJsonObject prv_obj, QString server)
{
    QString tmp_fname = CacheManager::globalInstance()->makeLocalFile(MLUtil::makeRandomId(10) + ".PrvManagerdlg.prv");
    TextFile::write(tmp_fname, QJsonDocument(prv_obj).toJson());
    CacheManager::globalInstance()->setTemporaryFileDuration(tmp_fname, 600);
    QString cmd = "prv";
    QStringList args;
    args << "ensure-remote" << tmp_fname << "--server=" + server;

    execute_command_in_separate_thread(cmd, args, q, SLOT(slot_restart_thread()));

    /*
    if (QProcess::execute(cmd, args) == 0) {
        restart_thread();
    }
    else {
        QMessageBox::warning(0,"Ensure local","Problem running: "+cmd+" "+args.join(" "));
    }
    */
}

void PrvManagerDialogPrivate::restart_thread()
{
    if (m_thread.isRunning()) {
        m_thread.requestInterruption();
        m_thread.wait(5000);
        m_thread.terminate();
    }
    m_thread.prv_objects = m_prv_objects;
    m_thread.server_names = m_server_names;
    m_thread.start();
}

QString exec_process_and_return_output(QString cmd, QStringList args)
{
    QProcess P;
    P.setReadChannelMode(QProcess::MergedChannels);
    P.start(cmd, args);
    P.waitForStarted();
    P.waitForFinished(-1);
    return P.readAll();
}

bool PrvManagerDialogThread::check_if_on_local_disk(QJsonObject prv_obj)
{
    QString cmd = "prv";
    QStringList args;
    args << "locate";
    args << "--checksum=" + prv_obj["original_checksum"].toString();
    args << "--fcs=" + prv_obj["original_fcs"].toString();
    args << QString("--size=%1").arg(prv_obj["original_size"].toVariant().toLongLong());
    args << QString("--original_path=%1").arg(prv_obj["original_path"].toVariant().toLongLong());
    args << "--local-only";
    QString output = exec_process_and_return_output(cmd, args);
    return !output.isEmpty();
}

bool PrvManagerDialogThread::check_if_on_server(QJsonObject prv_obj, QString server_name)
{
    QString cmd = "prv";
    QStringList args;
    args << "locate";
    args << "--checksum=" + prv_obj["original_checksum"].toString();
    args << "--fcs=" + prv_obj["original_fcs"].toString();
    args << QString("--size=%1").arg(prv_obj["original_size"].toVariant().toLongLong());
    args << "--server=" + server_name;
    QString output = exec_process_and_return_output(cmd, args);
    return !output.isEmpty();
}

void PrvManagerDialogThread::run()
{
    TaskProgress task("PRV Upload info...");
    QStringList names = prv_objects.keys();
    for (int i = 0; i < names.count(); i++) {
        task.setProgress((i + 0.5) / names.count());
        QString name = names[i];
        if (QThread::currentThread()->isInterruptionRequested())
            return;
        QJsonObject obj = prv_objects[name];
        {
            task.log() << "check if on local disk" << name;
            bool is = check_if_on_local_disk(obj);
            {
                QMutexLocker locker(&results_mutex);
                if (is)
                    results[name].on_local_disk = YES;
                else
                    results[name].on_local_disk = NO;
            }
            emit results_updated();
        }
        foreach (QString server_name, server_names) {
            task.log() << "check if on server" << name << server_name;
            if (QThread::currentThread()->isInterruptionRequested())
                return;
            bool is = check_if_on_server(obj, server_name);
            {
                QMutexLocker locker(&results_mutex);
                if (is)
                    results[name].on_server[server_name] = YES;
                else
                    results[name].on_server[server_name] = NO;
            }
            emit results_updated();
        }
    }
}
