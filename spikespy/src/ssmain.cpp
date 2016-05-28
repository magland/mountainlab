#include <QApplication>
#include <QScriptEngine>
#include <QDebug>
#include "sscontroller.h"
#include <QStringList>
#include "textfile.h"
#include "usagetracking.h"
#include "closemehandler.h"

/*
 * TO DO:
 * Clean up temporary files
 * */


int main(int argc, char *argv[]) {
	QApplication a(argc, argv);
	//MainWindow w;
	//w.show();
	CloseMeHandler::start();

	QStringList args;
	for (int i=1; i<argc; i++) {
		args << QString(argv[i]);
	}

	//args << "/tmp/tp903426a1_f13b_43a2_a1b8_5ff229f26018_run_java_script.js";

	bool connect_zeros=true;

	QString script_path;
	QList<QString> mda_paths;
	foreach (QString str,args) {
		if (str.endsWith(".js")) {
			script_path=str;
		}
		else if (str.endsWith(".mda")) {
			mda_paths << str;
		}
		else if (str=="do-not-connect-zeros") {
			connect_zeros=false;
		}
	}

	QScriptEngine *engine=new QScriptEngine;

	//qScriptRegisterMetaType(engine, myObjectToScriptValue, myObjectFromScriptValue);

	SSController SPIKESPY;
	QScriptValue SPIKESPY_value = engine->newQObject(&SPIKESPY);
	   //QScriptValue SPIKESPY_value = engine->newQObject(new QObject());
	engine->globalObject().setProperty("SPIKESPY", SPIKESPY_value);

	QString script;
	if (!script_path.isEmpty()) {
		script=read_text_file(script_path);
	} else {
		QList<QString> TL_paths;
		if (mda_paths.isEmpty()) {
			mda_paths << "/home/magland/tmp/ssview_spike_test_X.mda";
			mda_paths << "/home/magland/tmp/ssview_spike_test_X.mda";
			TL_paths << "/home/magland/tmp/ssview_spike_test_TL.mda";
			TL_paths << "/home/magland/tmp/ssview_spike_test_TL.mda";
			//mda_path="/tmp/1b781243ab416c14595f20d260a8c35a_run_java_script_variable_X.mda";
		}
		script="";
		script+="var W=SPIKESPY.createTimeSeriesWidget();\n\n";
		for (int i=0; i<mda_paths.count(); i++) {
			script+=QString("var V%1=SPIKESPY.createTimeSeriesView();\n").arg(i);
			script+=QString("var X%1=SPIKESPY.loadArray('%2');\n").arg(i).arg(mda_paths[i]);
			script+=QString("V%1.setData(X%1,true);\n").arg(i);
			if (!connect_zeros) {
				script+=QString("V%1.setConnectZeros(false);\n").arg(i);
			}
			script+=QString("V%1.setTitle(\"This is a test %1.\");\n").arg(i);
			script+=QString("W.addView(V%1);\n").arg(i);

			if (i<TL_paths.count()) {
				script+=QString("var TL%1=SPIKESPY.readArray('%2');\n").arg(i).arg(TL_paths[i]);
				script+=QString("V%1.setLabels(TL%1,true); V%1.initialize();\n").arg(i);
				//script+=QString("var L%1=SPIKESPY.createLabelView();\n").arg(i);
				//script+=QString("L%1.setLabels(TL%1,true);\n").arg(i);
				//script+=QString("W.addView(L%1);\n").arg(i);
			}

			script+="\n";
		}
	}

	QScriptValue result = engine->evaluate(script);
	if (result.isError()) {
		qWarning() << "Error running script: "+result.toString();
	}

    CleanupObject cleanup_object;
    QObject::connect(&a, SIGNAL(aboutToQuit()), &cleanup_object, SLOT(closing()));

	int ret=a.exec();

	engine->collectGarbage();

	delete engine;
	printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n",jnumfilesopen(),jmalloccount(),(int)jbytesallocated()*1.0/1000000);

	return ret;
}
