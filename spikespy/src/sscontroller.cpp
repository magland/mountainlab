#include "sscontroller.h"
#include "sstimeserieswidget_prev.h"
#include "sscommon.h"
#include "diskarraymodel_prev.h"
#include <QDateTime>
#include <QDebug>
#include <QDir>
#include <QFileInfo>
#include "sstimeseriesview_prev.h"
#include "sslabelview.h"
#include "mdaobject.h"
#include "diskreadmdaold.h"

SSController::SSController()
{

}

SSController::~SSController()
{

}

QWidget *SSController::createTimeSeriesWidget() {
	SSTimeSeriesWidget *W=new SSTimeSeriesWidget();
	W->setAttribute(Qt::WA_DeleteOnClose);
	W->showNormal();
	W->resize(1000,500);
	W->move(300,300);
	return W;
}

QWidget *SSController::createTimeSeriesView() {
	SSTimeSeriesView *V=new SSTimeSeriesView();
	return V;
}

QWidget *SSController::createLabelView()
{
	SSLabelView *V=new SSLabelView();
	return V;
}

QObject *SSController::loadArray(QString path) {
	SSARRAY *X=new SSARRAY();
	X->setPath(path.toLatin1().data());

	return X;
}

QObject *SSController::readArray(QString path)
{
    DiskReadMdaOld *ret=new DiskReadMdaOld;
	ret->setPath(path);
	return ret;
}

static QList<QString> s_paths_to_remove;

void CleanupObject::closing()
{
    for (int i=0; i<s_paths_to_remove.count(); i++) {
        QString path=s_paths_to_remove[i];

        QDateTime time=QFileInfo(path).lastModified();
        QString timestamp=time.toString("yyyy-mm-dd-hh-mm-ss");
        QString tmp=QFileInfo(path).path()+"/spikespy."+QFileInfo(path).completeBaseName()+"."+timestamp;
        if (QDir(tmp).exists()) {
            qDebug()  << "Removing directory: " << tmp;
            QStringList list=QDir(tmp).entryList(QStringList("*.mda"),QDir::Files|QDir::NoDotAndDotDot);
            foreach (QString A,list) {
                QFile::remove(tmp+"/"+A);
            }
            QDir(QFileInfo(tmp).path()).rmdir(QFileInfo(tmp).fileName());
        }

        qDebug()  << "Removing file: " << path;
        QFile(path).remove();
    }
}

void removeOnClose(const QString &path)
{
    s_paths_to_remove << path;
}
