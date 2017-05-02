/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/29/2017
*******************************************************/
#ifndef INITIALIZE_CONFUSION_MATRIX_H
#define INITIALIZE_CONFUSION_MATRIX_H

#include <QThread>
#include <mountainprocessrunner.h>

class Initialize_confusion_matrix : public QThread {
    Q_OBJECT
public:
    Initialize_confusion_matrix()
    {
        QObject::connect(this, SIGNAL(finished()), this, SIGNAL(finishedInGui())); //so we can use lambda expressions
    }

    // input
    bool relabel = false;
    QString firings1;
    QString firings2;

    // output
    QString confusion_matrix;
    QString matched_firings;
    QString label_map;
    QString firings2_relabeled;

    void run()
    {
        MountainProcessRunner MPR;
        MPR.setProcessorName("mountainsort.confusion_matrix");
        QVariantMap params;
        params["firings1"] = firings1;
        params["firings2"] = firings2;
        params["match_matching_offset"] = 30;
        if (relabel)
            params["relabel_firings2"] = "true";
        //params["_force_run"] = "true";
        MPR.setInputParameters(params);
        confusion_matrix = MPR.makeOutputFilePath("confusion_matrix_out");
        matched_firings = MPR.makeOutputFilePath("matched_firings_out");
        label_map = MPR.makeOutputFilePath("label_map_out");
        if (relabel)
            firings2_relabeled = MPR.makeOutputFilePath("firings2_relabeled_out");
        else
            firings2_relabeled = firings2;
        MPR.runProcess();
    }
signals:
    void finishedInGui(); //so we can use lambda expressions
};

#endif // INITIALIZE_CONFUSION_MATRIX_H
