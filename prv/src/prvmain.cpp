#include <QFile>
#include <QDebug>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTime>
#include <QCoreApplication>
#include <QThread>
#include <QDir>
#include <QJsonArray>
#include <QSettings>
#include <QProcess>
#include <QUrl>
//#include <QStandardPaths>
#include <QHostInfo>
#include <QCommandLineParser>
#include <QUrlQuery>
#include "cachemanager.h"
#include "prvfile.h"
#include "mlcommon.h"
#include "mlnetwork.h"

QString get_tmp_path();
QString get_server_url(QString url_or_server_name);
QStringList get_local_search_paths();
QJsonArray get_remote_servers(QString server_name = "");

void print(QString str);
void println(QString str);
bool is_file(QString path);
bool is_folder(QString path);

QString get_user_name()
{
#ifdef Q_OS_LINUX
    return qgetenv("USER");
#else
    return "user-name-not-implemented-yet-on-non-linux";
// WW: not a standard way to fetch username
//QStringList home_path = QStandardPaths::standardLocations(QStandardPaths::HomeLocation);
//return home_path.first().split(QDir::separator()).last();
#endif
}

namespace MLUtils {
class ApplicationCommand {
public:
    enum { ShowHelp = -1000 };
    virtual QString commandName() const = 0;
    virtual QString description() const { return QString(); }
    virtual void prepareParser(QCommandLineParser&) {}
    virtual ~ApplicationCommand() {}
    virtual int execute(const QCommandLineParser&) { return 0; }
};

class ApplicationCommandParser {
public:
    ApplicationCommandParser() {}
    bool process(const QCoreApplication& app)
    {
        m_result = 0;
        QCommandLineParser parser;
        parser.addPositionalArgument("command", "Command to execute");
        parser.addHelpOption();
        // TODO: build description from list of commands
        QString desc;
        const QLatin1Char nl('\n');
        desc += nl;
        desc += "Commands:" + nl;
        int longestCommandName = 0;
        foreach (ApplicationCommand* cmd, m_commands) {
            longestCommandName = qMax(longestCommandName, cmd->commandName().size());
        }

        foreach (ApplicationCommand* cmd, m_commands) {
            desc += cmd->commandName();
            if (cmd->description().isEmpty()) {
                desc += nl;
            }
            else {
                int spaces = longestCommandName - cmd->commandName().size() + 1;
                desc += QString(spaces, ' ');
                desc += '\t' + cmd->description() + nl;
            }
        }
        parser.setApplicationDescription(desc);
        parser.parse(app.arguments());
        if (parser.positionalArguments().isEmpty()) {
            parser.showHelp(1);
            return false;
        }
        QString command = parser.positionalArguments().first();
        foreach (ApplicationCommand* cmd, m_commands) {
            if (cmd->commandName() == command) {
                parser.clearPositionalArguments();
                parser.setApplicationDescription(QString());
                cmd->prepareParser(parser);
                parser.process(app);
                m_result = cmd->execute(parser);
                if (m_result == ApplicationCommand::ShowHelp)
                    parser.showHelp(m_result);
                return true;
            }
        }
        // command not found
        parser.showHelp(0);
        return false;
    }

    ~ApplicationCommandParser()
    {
        qDeleteAll(m_commands);
    }

    void addCommand(ApplicationCommand* cmd)
    {
        m_commands << cmd;
    }
    int result() const { return m_result; }

private:
    QList<ApplicationCommand*> m_commands;
    int m_result = 0;
};
}

namespace PrvCommands {

/*
class Sha1SumCommand : public MLUtils::ApplicationCommand {
public:
    QString commandName() const { return QStringLiteral("sha1sum"); }
    QString description() const { return QStringLiteral("Calculate sha1 hash for file"); }
    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("file_name", "File name");
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name
        if (args.isEmpty()) {
            return ShowHelp; // show help
        }
        QString path = args.first();
        if (!QFile::exists(path)) {
            qWarning() << "No such file: " + path;
            return -1;
        }
        return sha1sum(path, QVariantMap());
    }

private:
    int sha1sum(QString path, const QVariantMap& params) const
    {
        Q_UNUSED(params)
        QString checksum = MLUtil::computeSha1SumOfFile(path);
        if (checksum.isEmpty())
            return -1;
        println(checksum);
        return 0;
    }
};
*/

class StatCommand : public MLUtils::ApplicationCommand {
public:
    QString commandName() const { return "stat"; }
    QString description() const { return "Stat a file"; }

    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("file_name", "File name");
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name
        if (args.isEmpty()) {
            return ShowHelp; // show help
        }
        QString path = args.first();
        if (!QFile::exists(path)) {
            qWarning() << "No such file: " + path;
            return -1;
        }
        return stat(path, QVariantMap());
    }

private:
    int stat(QString path, const QVariantMap& params) const
    {
        Q_UNUSED(params)
        QString checksum = MLUtil::computeSha1SumOfFile(path);
        if (checksum.isEmpty()) {
            QJsonObject obj;
            obj["error"] = "checksum is empty for " + path;
            println(QJsonDocument(obj).toJson());
            return -1;
        }
        {
            QJsonObject obj;
            obj["checksum"] = checksum;
            obj["fcs"] = "head1000-" + MLUtil::computeSha1SumOfFileHead(path, 1000);
            obj["size"] = QFileInfo(path).size();
            println(QJsonDocument(obj).toJson());
            return 0;
        }
    }
};

/*
class ListSubserversCommand : public MLUtils::ApplicationCommand {
public:
    QString commandName() const { return QStringLiteral("list-subservers"); }
    QString description() const { return QStringLiteral("List sub servers"); }

    int execute(const QCommandLineParser&)
    {
        return list_subservers();
    }

private:
    int list_subservers() const
    {
        QJsonArray remote_servers = get_remote_servers();
        for (int i = 0; i < remote_servers.count(); i++) {
            QJsonObject server0 = remote_servers[i].toObject();
            QString host = server0["host"].toString();
            int port = server0["port"].toInt();
            QString url_path = server0["path"].toString();
            QString url0 = host + ":" + QString::number(port) + url_path + QString("/?a=list-subservers");
            url0 += "&passcode=" + server0["passcode"].toString();
            println("Connecting to " + url0);
            QString txt = http_get_text_curl_0(url0);
            print(txt + "\n\n");
        }
        return 0;
    }
};
*/

/*
class UploadCommand : public MLUtils::ApplicationCommand {
public:
    UploadCommand() {}
    QString commandName() const { return QStringLiteral("upload"); }
    QString description() const { return QStringLiteral("Upload files to server"); }
    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument(QStringLiteral("path"), QStringLiteral("File or directory name"));
        parser.addPositionalArgument(QStringLiteral("server"), QStringLiteral("Server name or URL"));
        parser.addOption(QCommandLineOption(
            QStringList() << QStringLiteral("param") << QStringLiteral("p"),
            QStringLiteral("Custom parameter which will be set while uploading"),
            QStringLiteral("param")));

        parser.addOption(QCommandLineOption(
            QStringList() << QStringLiteral("include") << QStringLiteral("i"),
            QStringLiteral("file patterns to include in upload (can be given multiple times)"), "include"));
        parser.addOption(QCommandLineOption(
            QStringList() << QStringLiteral("exclude") << QStringLiteral("e"),
            QStringLiteral("file patterns to exclude from upload (can be given multiple times)"), "exclude"));
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name
        if (args.size() < 2) {
            println("Too few arguments to prv upload.");
            qDebug().noquote() << args;
            return -1;
        }

        QString src_path = args.value(0);
        QString server_url = get_server_url(args.value(1));
        if (src_path.isEmpty()) {
            return ShowHelp;
        }
        if (!QFile::exists(src_path)) {
            qWarning() << "No such file: " + src_path;
            return -1;
        }
        if (src_path.endsWith(".prv")) {
            qWarning() << "Cannot upload a .prv file.";
            return -1;
        }
        // parser custom options
        QVariantMap customParams;
        if (parser.isSet("param")) {
            foreach (const QString param, parser.values("param")) {
                int idx = param.indexOf('=');
                if (idx <= 0)
                    continue; // TODO: error?
                QString key = param.left(idx);
                QVariant value = param.mid(idx + 1);
                customParams[key] = value;
            }
        }
        QStringList includes = parser.values("include");
        QStringList excludes = parser.values("exclude");
        m_totalSize = 0;
        upload(src_path, server_url, customParams, includes, excludes);
        return 0;
    }

private:
    int64_t upload(const QString& path, const QString& server_url, const QVariantMap& params,
        const QStringList& includes = QStringList(), const QStringList& excludes = QStringList())
    {
        int64_t totalSize = 0;
        QFileInfo pathInfo = path;
        if (!pathInfo.isDir()) {
            return uploadFile(path, server_url, params);
        }
        QDir dir = path;
        QFileInfoList entries = dir.entryInfoList(includes, QDir::Files, QDir::Name);
        foreach (const QFileInfo& fileInfo, entries) {
            bool exclude = false;
            foreach (const QString& exPattern, excludes) {
                QRegExp rx(exPattern);
                rx.setPatternSyntax(QRegExp::WildcardUnix);
                if (rx.exactMatch(fileInfo.fileName())) {
                    exclude = true;
                }
                if (exclude)
                    break;
            }
            if (exclude)
                continue;
            int64_t res = uploadFile(fileInfo.absoluteFilePath(), server_url, params);
            if (res > 0) {
                totalSize += res;
            }
        }
        // recursive into subdirectories unless explicitly excluded
        entries = dir.entryInfoList(QStringList(), QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
        foreach (const QFileInfo& fileInfo, entries) {
            bool exclude = false;
            foreach (const QString& exPattern, excludes) {
                QRegExp rx(exPattern);
                rx.setPatternSyntax(QRegExp::WildcardUnix);
                if (rx.exactMatch(fileInfo.fileName())) {
                    exclude = true;
                }
                if (exclude)
                    break;
            }
            if (exclude)
                continue;
            int64_t res = upload(fileInfo.absoluteFilePath(), server_url, params, includes, excludes);
            if (res > 0) {
                totalSize += res;
            }
        }
        return totalSize;
    }

    int64_t uploadFile(const QString& path, const QString& server_url, const QVariantMap& params)
    {
        QFileInfo finfo = path;
        int size0 = finfo.size();
        if (size0 == 0) {
            println("File is empty... skipping: " + path);
            return 0;
        }
        QString checksum00 = MLUtil::computeSha1SumOfFile(path);
        if (checksum00.isEmpty()) {
            println("checksum is empty for file: " + path);
            return -1;
        }
        QUrl url(server_url);
        QUrlQuery query;
        //query.addQueryItem("a", "upload");
        query.addQueryItem("checksum", checksum00);
        query.addQueryItem("size", QString::number(size0));
        QJsonObject info;
        info["src_path"] = QDir::current().absoluteFilePath(path);
        info["server_url"] = server_url;
        info["user"] = get_user_name();
        info["local_host_name"] = QHostInfo::localHostName();
        info["date_uploaded"] = QDateTime::currentDateTime().toString("yyyy-MM-dd:hh-mm-ss");
        info["params"] = QJsonObject::fromVariantMap(params);
        QString info_json = QJsonDocument(info).toJson();
        query.addQueryItem("info", QUrl::toPercentEncoding(info_json.toUtf8()));
        url.setQuery(query);

        //QString ret=MLNetwork::httpPostFileSync(path,url.toString());
        QString ret = MLNetwork::httpPostFileParallelSync(path, url.toString());

//        QString ret = NetUtils::httpPostFile(url, path, [](qint64 bytesSent, qint64 bytesTotal) {
//                if (bytesTotal == 0) {
//                printf("\33[2K\r");
//                return;
//                }
//                printf("\33[2K\r");
//                unsigned int prc = qRound(50*((double)bytesSent/(double)bytesTotal));
//                for(unsigned int i = 0; i < prc; ++i) printf("#");
//                for(unsigned int i = prc; i < 50; ++i) printf(".");
//                fflush(stdout);
//        });

        if (ret.isEmpty()) {
            qWarning() << "Problem posting file to: " + url.toString();
            return -1;
        }
        QJsonParseError err0;
        QJsonObject obj = QJsonDocument::fromJson(ret.toUtf8(), &err0).object();
        if (err0.error != QJsonParseError::NoError) {
            println(QString("Error uploading file."));
            return -1;
        }
        if (!obj["success"].toBool()) {
            QString error = obj["error"].toString();
            println(QString("Error uploading file: %1").arg(error));
            return -1;
        }
        m_totalSize += size0;
        println(QString("Uploaded file %1 (%2 MB, total %3 MB).").arg(finfo.fileName()).arg(size0 * 1.0 / 1e6).arg(m_totalSize * 1.0 / 1e6));
        return size0;
    }
    int64_t m_totalSize = 0;
};
*/

class CreateCommand : public MLUtils::ApplicationCommand {
public:
    QString commandName() const { return "create"; }
    QString description() const { return "Creates a new prv file"; }
    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("source", "Source file or directory name");
        parser.addPositionalArgument("dest", "Destination file or directory name", "[dest]");
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name
        if (args.isEmpty()) {
            return ShowHelp; // show help
        }

        QString src_path = args.value(0);
        QString dst_path = args.value(1);
        if (src_path.isEmpty()) {
            return -1;
        }
        if (!QFile::exists(src_path)) {
            qWarning() << "No such file: " + src_path;
            return -1;
        }
        if (dst_path.isEmpty()) {
            dst_path = QFileInfo(src_path).fileName() + ".prv";
        }
        if (!dst_path.endsWith(".prv")) {
            println("Destination file must end with .prv");
            return -1;
        }
        QVariantMap params;
        if (is_file(src_path)) {
            int ret = create_file_prv(src_path, dst_path, params);
            if (ret != 0)
                return ret;
            return 0;
        }
        else if (is_folder(src_path)) {
            return create_folder_prv(src_path, dst_path, params);
        }
        else {
            qWarning() << "not sure why file is not a file nor a folder.";
            return -1;
        }
    }

private:
    int create_file_prv(QString src_path, QString dst_path, const QVariantMap& params) const
    {
        println("making prv file: " + dst_path);
        PrvFile PF;
        PrvFileCreateOptions opts;
        opts.create_temporary_files = params.contains("create-temporary-files");
        PF.createFromFile(src_path, opts);
        if (!PF.write(dst_path))
            return -1;
        return 0;
    }

    int create_folder_prv(QString src_path, QString dst_path, const QVariantMap& params) const
    {
        println("making folder prv: " + src_path);
        PrvFile PF;
        PrvFileCreateOptions opts;
        opts.create_temporary_files = params.contains("create-temporary-files");
        PF.createFromDirectory(src_path, opts);
        if (!PF.write(dst_path))
            return -1;
        return 0;
    }
};

/*
class RecoverCommand : public MLUtils::ApplicationCommand {
public:
    QString commandName() const { return "recover"; }

    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("source", "Source PRV file name");
        parser.addPositionalArgument("dest", "Destination file or directory name", "[dest]");
        parser.addOption(QCommandLineOption("recover-all-prv-files", ""));
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name
        if (args.isEmpty()) {
            return ShowHelp; // show help
        }
        QString src_path = args.value(0);
        QString dst_path = args.value(1);
        if (src_path.isEmpty()) {
            return ShowHelp;
        }
        if (!QFile::exists(src_path)) {
            qWarning() << "No such file: " + src_path;
            return -1;
        }
        if (!src_path.endsWith(".prv")) {
            println("prv file must have .prv extension");
            return -1;
        }
        if (dst_path.isEmpty()) {
            QString f0 = QFileInfo(src_path).fileName();
            dst_path = f0.left(f0.count() - 4); //remove .prv extension
        }
        PrvFile prv_file(src_path);
        PrvFileRecoverOptions opts;
        opts.recover_all_prv_files = parser.isSet("recover-all-prv-files");
        opts.locate_opts.local_search_paths = get_local_search_paths();
        opts.locate_opts.search_remotely = true;
        opts.locate_opts.remote_servers = get_remote_servers();
        if (prv_file.representsFile()) {
            if (!prv_file.recoverFile(dst_path, opts))
                return -1;
        }
        else {
            if (!prv_file.recoverFolder(dst_path, opts))
                return -1;
        }
        return 0;
    }
};
*/

class LocateOrDownloadCommand : public MLUtils::ApplicationCommand {
public:
    LocateOrDownloadCommand(const QString& cmd)
        : m_cmd(cmd)
    {
    }

public:
    QString commandName() const { return m_cmd; }
    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("file_name", "PRV file name", "[file_name]");
        if (m_cmd == "download") {
            parser.addPositionalArgument("output_file_name", "output file name", "[output_file_name]");
        }
        // --checksum=[] --fcs=[optional] --size=[]
        parser.addOption(QCommandLineOption("checksum", "checksum", "[]"));
        parser.addOption(QCommandLineOption("fcs", "fcs", "[optional]"));
        parser.addOption(QCommandLineOption("original_path", "original_path", "[optional]"));
        parser.addOption(QCommandLineOption("size", "size", "[]"));
        parser.addOption(QCommandLineOption("search_path", "search_path", "path to search (leave empty to search default locations)"));
        parser.addOption(QCommandLineOption("server", "server", "[url or name of prvbucket server to search]"));
        parser.addOption(QCommandLineOption("verbose", "verbose"));
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name

        bool verbose = parser.isSet("verbose");

        QJsonObject obj;
        if (parser.isSet("checksum")) {
            obj["original_checksum"] = parser.value("checksum");
            obj["original_fcs"] = parser.value("fcs");
            obj["original_size"] = parser.value("size").toLongLong();
            obj["original_path"] = parser.value("original_path");
        }
        else {
            QString src_path = args.value(0);
            if (src_path.isEmpty()) {
                println("Source path is empty");
                return -1;
            }
            if (!QFile::exists(src_path)) {
                qWarning() << "No such file: " + src_path;
                return -1;
            }
            if (src_path.endsWith(".prv")) {
                obj = QJsonDocument::fromJson(TextFile::read(src_path).toUtf8()).object();
            }
            else {
                if (m_cmd == "locate") {
                    obj["original_checksum"] = MLUtil::computeSha1SumOfFile(src_path);
                    obj["original_fcs"] = "head1000-" + MLUtil::computeSha1SumOfFileHead(src_path, 1000);
                    obj["original_size"] = QFileInfo(src_path).size();
                }
                else {
                    println("prv file must have .prv extension");
                    return -1;
                }
            }
        }
        if (obj.contains("original_checksum")) {
            QVariantMap params;
            if (parser.isSet("search_path"))
                params.insert("search_path", parser.value("search_path"));
            if (parser.isSet("server"))
                params.insert("server", parser.value("server"));
            QString fname_or_url = locate_file(obj, params, verbose);
            if (fname_or_url.isEmpty()) {
                println("Unable to locate file");
                return -1;
            }
            if (m_cmd == "locate") {
                println(fname_or_url);
            }
            else {
                println("download: " + fname_or_url);
                QString dst_fname = args.value(1);
                QString cmd;
                if (is_url(fname_or_url)) {
                    QString tmp = "";
                    if (!dst_fname.isEmpty()) {
                        if (QFileInfo(dst_fname).isDir())
                            tmp += "-P " + dst_fname;
                        else
                            tmp += "-O " + dst_fname;
                    }
                    cmd = QString("wget %1 %2").arg(fname_or_url).arg(tmp);
                }
                else {
                    QString tmp = "";
                    if (!dst_fname.isEmpty())
                        tmp += "> " + dst_fname;
                    cmd = QString("cat %1 %2").arg(fname_or_url).arg(tmp);
                }
                println(QString("Running: %1").arg(cmd));
                return system(cmd.toUtf8().data());
            }
        }
        else {
            //it must be a directory
            locate_directory(obj, verbose);
        }
        return 0;
    }

private:
    QString m_cmd;

    QString locate_file(const QJsonObject& obj, const QVariantMap& params, bool verbose) const
    {
        PrvFile prvf(obj);
        PrvFileLocateOptions opts;
        if (params.contains("server")) {
            opts.search_locally = false;
            opts.search_remotely = true;
            opts.remote_servers << params["server"].toString();
        }
        else if (params.contains("search_path")) {
            opts.search_locally = true;
            opts.search_remotely = false;
            opts.local_search_paths.clear();
            opts.local_search_paths << params["search_path"].toString();
        }
        else {
            opts.local_search_paths = get_local_search_paths();
            opts.search_remotely = false;
        }

        opts.verbose = verbose;
        QString fname_or_url = prvf.locate(opts);
        if (fname_or_url.isEmpty())
            return "";
        return fname_or_url;
    }

    int locate_directory(const QJsonObject& obj, bool verbose) const
    {
        PrvFile prvf(obj);
        PrvFileLocateOptions opts;
        opts.search_locally = true;
        opts.search_remotely = false;
        opts.verbose = verbose;
        opts.local_search_paths = get_local_search_paths();
        QString fname_or_url = prvf.locateDirectory(opts);
        if (fname_or_url.isEmpty())
            return -1;
        println(fname_or_url);
        return 0;
    }

    /*
    int download_file(const QJsonObject& obj, const QVariantMap& params) const
    {
        PrvFile prvf(obj);
        PrvFileLocateOptions opts;
        opts.local_search_paths = get_local_search_paths();
        opts.search_remotely = true;
        opts.remote_servers = get_remote_servers();

        if (params.contains("path")) {
            opts.local_search_paths.clear();
            opts.local_search_paths << params["path"].toString();
            opts.search_remotely = false;
        }

        QString fname_or_url = prvf.locate(opts);
        if (fname_or_url.isEmpty())
            return -1;

        QString cmd;
        if (is_url(fname_or_url)) {
            cmd = QString("curl %1").arg(fname_or_url);
        }
        else {
            cmd = QString("cat %1").arg(fname_or_url);
        }
        return system(cmd.toUtf8().data());
    }
    */
};

/*
QList<QJsonObject> get_all_input_prv_objects(QJsonObject obj)
{
    QList<QJsonObject> ret;
    ret << obj;
    QJsonArray processes = obj["processes"].toArray();
    for (int i = 0; i < processes.count(); i++) {
        QJsonObject process0 = processes[i].toObject();
        QJsonObject inputs = process0["inputs"].toObject();
        QStringList names = inputs.keys();
        foreach (QString name, names) {
            ret << get_all_input_prv_objects(inputs[name].toObject());
        }
    }
    return ret;
}

QList<QJsonObject> get_all_output_prv_objects(QJsonObject obj)
{
    QList<QJsonObject> ret;
    QJsonArray processes = obj["processes"].toArray();
    for (int i = 0; i < processes.count(); i++) {
        QJsonObject process0 = processes[i].toObject();
        QJsonObject outputs = process0["outputs"].toObject();
        QStringList names = outputs.keys();
        foreach (QString name, names) {
            ret << outputs[name].toObject();
            ret << get_all_output_prv_objects(outputs[name].toObject());
        }
    }
    return ret;
}

QList<QJsonObject> get_all_inputs_that_are_not_outputs(const QList<QJsonObject>& inputs, const QList<QJsonObject>& outputs)
{
    QList<QJsonObject> ret;

    QSet<QString> output_path_set;
    for (int i = 0; i < outputs.count(); i++) {
        output_path_set.insert(outputs[i]["original_path"].toString());
    }
    QSet<QString> used_input_path_set;
    for (int i = 0; i < inputs.count(); i++) {
        QString tmp = inputs[i]["original_path"].toString();
        if (!used_input_path_set.contains(tmp)) {
            if (!output_path_set.contains(tmp)) {
                ret << inputs[i];
                used_input_path_set.insert(tmp);
            }
        }
    }
    return ret;
}

QList<PrvFile> extract_raw_only_prv_files(const PrvFile& prv_file)
{
    QList<QJsonObject> input_objects = get_all_input_prv_objects(prv_file.object());
    QList<QJsonObject> output_objects = get_all_output_prv_objects(prv_file.object());
    QList<QJsonObject> raw_objects = get_all_inputs_that_are_not_outputs(input_objects, output_objects);

    if (raw_objects.count() == 1) {
        if (raw_objects[0]["original_checksum"].toString() == prv_file.checksum()) {
            QList<PrvFile> ret0;
            ret0 << prv_file;
            return ret0;
        }
    }

    QList<PrvFile> ret;
    foreach (QJsonObject raw_object, raw_objects) {
        ret << PrvFile(raw_object);
    }
    return ret;
}
*/

/*
class EnsureLocalRemoteCommand : public MLUtils::ApplicationCommand {
public:
    EnsureLocalRemoteCommand(QString cmd)
        : m_cmd(cmd)
    {
    }

public:
    QString commandName() const { return m_cmd; }
    QString description() const { return QStringLiteral("If the underlying file is not in the search paths of the local machine, copy it, download it, or regenerate it."); }
    void prepareParser(QCommandLineParser& parser)
    {
        parser.addPositionalArgument("file_name", "PRV file name", "[file_name]");
        //parser.addOption(QCommandLineOption("regenerate-if-needed", "recreate file using processing if needed"));
        parser.addOption(QCommandLineOption("download-if-needed", "download file from a remote server if needed"));
        //parser.addOption(QCommandLineOption("download-and-regenerate-if-needed", "download file from a remote server if needed"));
        parser.addOption(QCommandLineOption("raw-only", "if processing provenance is available, only upload the raw files."));
        parser.addOption(QCommandLineOption("server", "name or url of server", "[name|url]"));
    }
    int execute(const QCommandLineParser& parser)
    {
        QStringList args = parser.positionalArguments();
        args.removeFirst(); // remove command name

        PrvFile prv_file0(args.value(0));
        QList<PrvFile> prv_files;
        if (parser.isSet("raw-only")) {
            prv_files = extract_raw_only_prv_files(prv_file0);
        }
        else {
            prv_files << prv_file0;
        }
        foreach (PrvFile prv_file, prv_files) {
            println("Handling file " + prv_file.prvFilePath() + ": original_path=" + prv_file.originalPath());
            QString dst_path = CacheManager::globalInstance()->makeLocalFile(prv_file.checksum(), CacheManager::LongTerm);

            QString server; //for ensure-remote only
            if (m_cmd == "ensure-remote") {
                // Let's check whether it is on the server
                server = parser.value("server");
                if (server.isEmpty()) {
                    println("Missing required parameter: server");
                    return -1;
                }
                PrvFileLocateOptions opts;
                opts.search_locally = false;
                opts.search_remotely = true;
                opts.remote_servers = get_remote_servers(server);
                QString url = prv_file.locate(opts);
                if (!url.isEmpty()) {
                    println("File is already on server: " + url);
                    break;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////
            // Check whether it is already on the local system in a directory that
            // is searched. For example, the temporary directories are usually searched.
            PrvFileLocateOptions opts;
            opts.local_search_paths = get_local_search_paths();
            opts.search_remotely = false;
            QString path0 = prv_file.locate(opts);
            if (!path0.isEmpty()) {
                // Terrific!!
                if (m_cmd == "ensure-local") {
                    println("File is already on the local system: " + path0);
                    break;
                }
                else if (m_cmd == "ensure-remote") {
                    println("File is on the local system: " + path0);
                    if (upload_file_to_server(path0, server))
                        break;
                    else
                        return -1;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////////////
            // Next we check whether the file is in its original location
            // This is the one and only time we can actually look at the original path
            // That's because it's likely that this .prv file was just created a moment ago
            // If so, we simply can copy the file into the temporary directory, to ensure it
            // is on the local machine.
            QString path1 = prv_file.originalPath();
            if (QFile::exists(path1)) {
                if (MLUtil::computeSha1SumOfFile(path1) == prv_file.checksum()) {
                    if (QFile::exists(dst_path))
                        QFile::remove(dst_path);
                    if (QFile::copy(path1, dst_path)) {
                        println("Copied file to: " + dst_path);
                        if (m_cmd == "ensure-remote") {
                            if (upload_file_to_server(dst_path, server))
                                break;
                            else
                                return -1;
                        }
                        else
                            break;
                    }
                    else {
                        println("PROBLEM: Unable to copy file to: " + dst_path);
                        return -1;
                    }
                }
            }
            ////////////////////////////////////////////////////////////////////////////////

//            ////////////////////////////////////////////////////////////////////////////////
//            // Next we see whether we can recreate it using local processing, if this option
//            // has been enabled
//            if (parser.isSet("regenerate-if-needed")) {
//                println("Attempting to regenerate file.");
//                QString tmp_path = resolve_prv_file(prv_file.prvFilePath(), false, true);
//                if (!tmp_path.isEmpty()) {
//                    println("Generated file: " + tmp_path);
//                    if (m_cmd == "ensure-remote") {
//                        if (upload_file_to_server(dst_path, server))
//                            break;
//                        else
//                            return -1;
//                    }
//                    else
//                        break;
//                }
//            }
//            ////////////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////////////
            // Next we see whether we can download it if this option has been enabled
            if ((parser.isSet("download-if-needed")) || (!parser.value("server").isEmpty())) {
                println("Attempting to download file.");
                PrvFileRecoverOptions opts;
                opts.locate_opts.remote_servers = get_remote_servers();
                if (!parser.value("server").isEmpty()) {
                    opts.locate_opts.remote_servers = get_remote_servers(parser.value("server"));
                }
                opts.locate_opts.search_locally = false;
                opts.locate_opts.search_remotely = true;
                opts.recover_all_prv_files = false;
                if (prv_file.recoverFile(dst_path, opts)) {
                    if (m_cmd == "ensure-remote") {
                        if (upload_file_to_server(dst_path, server))
                            break;
                        else
                            return -1;
                    }
                    else
                        break;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////

//            ////////////////////////////////////////////////////////////////////////////////
//            // Next we see whether we can use a combination of downloads and processing
//            // if this option has been enabled
//            if (parser.isSet("download-and-regenerate-if-needed")) {
//                println("Attempting to download and regenerate file.");
//                QString tmp_path = resolve_prv_file(prv_file.prvFilePath(), true, true);
//                if (!tmp_path.isEmpty()) {
//                    println("Generated file: " + tmp_path);
//                    if (m_cmd == "ensure-remote") {
//                        if (upload_file_to_server(dst_path, server))
//                            break;
//                        else
//                            return -1;
//                    }
//                    else
//                        break;
//                }
//            }
//            ////////////////////////////////////////////////////////////////////////////////

            println("Unable to find file on local machine. Original path: " + prv_file.originalPath());
            return -1;
        }

        return 0;
    }

private:
    QString m_cmd;

    bool upload_file_to_server(QString path, QString server)
    {
        QString cmd = "prv";
        QStringList args;
        args << "upload" << path << server;
        return (QProcess::execute(cmd, args) == 0);
    }
};
*/

} // namespace PrvCommands

void myMessageOutput(QtMsgType type, const QMessageLogContext& context, const QString& msg)
{
    Q_UNUSED(type)
    Q_UNUSED(context)
    Q_UNUSED(msg)
    // do not display anything!
    /*
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        fprintf(stderr, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtInfoMsg:
        fprintf(stderr, "Info: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtWarningMsg:
        fprintf(stderr, "Warning: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtCriticalMsg:
        fprintf(stderr, "Critical: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtFatalMsg:
        fprintf(stderr, "Fatal: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        abort();
    }
    */
}

bool contains_argv(int argc, char* argv[], QString str)
{
    for (int i = 0; i < argc; i++) {
        if (QString(argv[i]) == str)
            return true;
    }
    return false;
}

int main(int argc, char* argv[])
{

    if (!contains_argv(argc, argv, "--verbose"))
        qInstallMessageHandler(myMessageOutput);

    QCoreApplication app(argc, argv);

    CacheManager::globalInstance()->setLocalBasePath(get_tmp_path());

    MLUtils::ApplicationCommandParser cmdParser;
    //cmdParser.addCommand(new PrvCommands::Sha1SumCommand);
    cmdParser.addCommand(new PrvCommands::StatCommand);
    cmdParser.addCommand(new PrvCommands::CreateCommand);
    cmdParser.addCommand(new PrvCommands::LocateOrDownloadCommand("locate"));
    cmdParser.addCommand(new PrvCommands::LocateOrDownloadCommand("download"));
    //cmdParser.addCommand(new PrvCommands::RecoverCommand);
    //cmdParser.addCommand(new PrvCommands::ListSubserversCommand);
    //cmdParser.addCommand(new PrvCommands::UploadCommand);
    //cmdParser.addCommand(new PrvCommands::EnsureLocalRemoteCommand("ensure-local"));
    //cmdParser.addCommand(new PrvCommands::EnsureLocalRemoteCommand("ensure-remote"));
    if (!cmdParser.process(app)) {
        return cmdParser.result();
    }
    return cmdParser.result();
}

void print(QString str)
{
    printf("%s", str.toUtf8().constData());
}

void println(QString str)
{
    printf("%s\n", str.toUtf8().constData());
}

QString get_tmp_path()
{
    QString temporary_path = MLUtil::tempPath();
    if (temporary_path.isEmpty())
        return "";
    QDir(temporary_path).mkdir("prv");
    return temporary_path + "/prv";
}

QString make_temporary_file()
{
    QString file_name = MLUtil::makeRandomId(10) + ".tmp";
    return get_tmp_path() + "/" + file_name;
}

bool is_file(QString path)
{
    return QFileInfo(path).isFile();
}
bool is_folder(QString path)
{
    return QFileInfo(path).isDir();
}

QString get_server_url(QString url_or_server_name)
{
    QJsonArray remote_servers = MLUtil::configValue("prv", "servers").toArray();
    for (int i = 0; i < remote_servers.count(); i++) {
        QJsonObject server0 = remote_servers[i].toObject();
        if (server0["name"].toString() == url_or_server_name) {
            QString host = server0["host"].toString();
            int port = server0["port"].toInt();
            QString url_path = server0["path"].toString();
            QString url0 = host + ":" + QString::number(port) + url_path;
            return url0;
        }
    }
    return url_or_server_name;
}

QJsonArray get_remote_servers(QString server_name)
{
    QJsonArray remote_servers = MLUtil::configValue("prv", "servers").toArray();
    if (!server_name.isEmpty()) {
        QJsonArray ret;
        foreach (QJsonValue s, remote_servers) {
            if (s.toObject()["name"].toString() == server_name)
                ret << s;
        }
        return ret;
    }
    return remote_servers;
}
