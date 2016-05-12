/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MSPROCESSOR_H
#define MSPROCESSOR_H

#include <QString>
#include <QVariant>
#include <QDebug>

class MSProcessorPrivate;
class MSProcessor
{
public:
	friend class MSProcessorPrivate;
	MSProcessor();
	virtual ~MSProcessor();

	QString name();
	QString version();
    QString description();
	QStringList inputFileParameters() const;
	QStringList outputFileParameters() const;
	QStringList requiredParameters() const;
	QStringList optionalParameters() const;

	virtual bool check(const QMap<QString,QVariant> &params)=0;
	virtual bool run(const QMap<QString,QVariant> &params)=0;

protected:
	void setName(const QString &name);
	void setVersion(const QString &version);
    void setDescription(const QString &description);
	void setInputFileParameters(const QString &p1,const QString &p2="",const QString &p3="",const QString &p4="");
	void setOutputFileParameters(const QString &p1,const QString &p2="",const QString &p3="",const QString &p4="");
	void setRequiredParameters(const QString &p1,const QString &p2="",const QString &p3="",const QString &p4="");
	void setOptionalParameters(const QString &p1,const QString &p2="",const QString &p3="",const QString &p4="");
	bool checkParameters(const QMap<QString,QVariant> &params);
private:
	MSProcessorPrivate *d;
};

#endif // MSPROCESSOR_H
