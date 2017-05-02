/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 9/29/2016
*******************************************************/

#include <QApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <icounter.h>
#include <qprocessmanager.h>
#include "mlcommon.h"
#include "prvguimainwindow.h"
#include "signal.h"
#include "objectregistry.h"

void sig_handler(int signum)
{
    (void)signum;
    QProcessManager* manager = ObjectRegistry::getObject<QProcessManager>();
    if (manager) {
        manager->closeAll();
    }
    abort();
}

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    QApplication::setApplicationName("prv-gui");
    QApplication::setApplicationVersion("0.1");

    ObjectRegistry registry;
    CounterManager* counterManager = new CounterManager;
    registry.addAutoReleasedObject(counterManager);

    //The process manager
    QProcessManager* processManager = new QProcessManager;
    registry.addAutoReleasedObject(processManager);
    signal(SIGINT, sig_handler);
    signal(SIGKILL, sig_handler);
    signal(SIGTERM, sig_handler);

    QCommandLineParser parser;
    parser.setApplicationDescription("prv-gui");

    parser.process(app);
    QStringList args = parser.positionalArguments();
    QString fname = args.value(0);
    if (fname.isEmpty()) {
        printf("No file argument specified.\n");
        return -1;
    }

    QJsonArray servers0 = MLUtil::configValue("prv", "servers").toArray();
    QStringList server_names;
    foreach (QJsonValue server, servers0) {
        server_names << server.toObject()["name"].toString();
    }

    PrvGuiMainWindow W;
    if (!W.loadPrv(fname)) {
        return -1;
    }
    W.setServerNames(server_names);
    W.showMaximized();
    W.refresh();
    W.startAllSearches();

    return app.exec();
}
