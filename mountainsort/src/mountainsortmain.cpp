/******************************************************
**
** Copyright (C) 2016 by Jeremy Magland
**
** This file is part of the MountainSort C++ project
**
** Some rights reserved.
** See accompanying LICENSE and README files.
**
*******************************************************/

#include <QCoreApplication>
#include <stdio.h>
#include "commandlineparams.h"
#include "msprocessmanager.h"
#include <QDebug>
#include <QJsonDocument>
#include <QJsonObject>
#include "diskreadmda.h"
#include "unit_tests.h"
#include "textfile.h"
#include "cachemanager.h"
#include "mlutils.h"

void print_usage();
void list_processors(const MSProcessManager* PM);
bool run_process(MSProcessManager* PM, QJsonObject process);

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);

    setbuf(stdout, NULL);

    CLParams CLP = commandlineparams(argc, argv);

    MSProcessManager* PM = MSProcessManager::globalInstance();
    PM->loadDefaultProcessors();

    QString arg1 = CLP.unnamed_parameters.value(0);
    QString arg2 = CLP.unnamed_parameters.value(1);

    if (arg1 == "process") {
        if (arg2.isEmpty()) {
            print_usage();
            return -1;
        }
        QString json = read_text_file(arg2);
        if (json.isEmpty()) {
            printf("Unable to open file or file is empty: %s\n", arg2.toLatin1().data());
        }
        QJsonObject obj = QJsonDocument::fromJson(json.toLatin1()).object();
        if (!run_process(PM, obj))
            return -1;
        else
            return 0;
    } else if (arg1 == "list-processors") {
        list_processors(PM);
        return 0;
    } else if (arg1 == "detail-processors") {
        PM->printDetails();
        return 0;
    } else if (arg1 == "spec") {
        PM->printJsonSpec();
        return 0;
    }

    if (CLP.unnamed_parameters.isEmpty()) {
        print_usage();
        return -1;
    }

    if (CLP.unnamed_parameters.count() == 1) {
        QString processor_name = CLP.unnamed_parameters[0];
        QJsonObject process;
        process["processor_name"] = processor_name;
        QJsonObject parameters;
        QStringList keys = CLP.named_parameters.keys();
        foreach(QString key, keys)
        {
            parameters[key] = CLP.named_parameters[key].toString();
        }
        process["parameters"] = parameters;
        if (run_process(PM, process))
            return 0;
        else
            return -1;
    } else {
        printf("Unexpected number of unnamed parameters: %d\n", CLP.unnamed_parameters.count());
    }

    return 0;
}

void print_usage()
{
    printf("mountainsort process [process.json]\n");
    printf("mountainsort pipeline [pipeline.json]\n");
    printf("mountainsort list-processors\n");
    printf("mountainsort detail-processors\n");
    printf("mountainsort [processor_name] --[param1]=[value1] --[param2]=[value2] ...\n");
}

void list_processors(const MSProcessManager* PM)
{
    QString str = PM->usageString();
    printf("%s\n", str.toLatin1().data());
}

bool run_process(MSProcessManager* PM, QJsonObject process)
{
    QString processor_name = process["processor_name"].toString();
    QJsonObject parameters = process["parameters"].toObject();
    QStringList keys = parameters.keys();
    QMap<QString, QVariant> params;
    foreach(QString key, keys)
    {
        params[key] = parameters[key].toString();
    }

    if (!PM->checkAndRunProcess(processor_name, params)) {
        return false;
    }

    return true;
}
