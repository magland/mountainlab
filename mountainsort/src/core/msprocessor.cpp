#include "msprocessor.h"

class MSProcessorPrivate
{
public:
	MSProcessor *q;
	QString m_name;
	QString m_version;
    QString m_description;
	QStringList m_input_file_parameters;
	QStringList m_output_file_parameters;
	QStringList m_required_parameters;
	QStringList m_optional_parameters;
};

MSProcessor::MSProcessor() {
	d=new MSProcessorPrivate;
	d->q=this;
}

MSProcessor::~MSProcessor() {
	delete d;
}

QString MSProcessor::name()
{
	return d->m_name;
}

QString MSProcessor::version()
{
    return d->m_version;
}

QString MSProcessor::description()
{
    return d->m_description;
}

QStringList MSProcessor::inputFileParameters() const
{
	return d->m_input_file_parameters;
}

QStringList MSProcessor::outputFileParameters() const
{
	return d->m_output_file_parameters;
}

QStringList MSProcessor::requiredParameters() const
{
	return d->m_required_parameters;
}

QStringList MSProcessor::optionalParameters() const
{
	return d->m_optional_parameters;
}

void MSProcessor::setName(const QString &name)
{
	d->m_name=name;
}

void MSProcessor::setVersion(const QString &version)
{
    d->m_version=version;
}

void MSProcessor::setDescription(const QString &description)
{
    d->m_description=description;
}

void MSProcessor::setInputFileParameters(const QString &p1, const QString &p2, const QString &p3, const QString &p4)
{
	if (!p1.isEmpty()) d->m_input_file_parameters << p1;
	if (!p2.isEmpty()) d->m_input_file_parameters << p2;
	if (!p3.isEmpty()) d->m_input_file_parameters << p3;
	if (!p4.isEmpty()) d->m_input_file_parameters << p4;
}

void MSProcessor::setOutputFileParameters(const QString &p1, const QString &p2, const QString &p3, const QString &p4)
{
	if (!p1.isEmpty()) d->m_output_file_parameters << p1;
	if (!p2.isEmpty()) d->m_output_file_parameters << p2;
	if (!p3.isEmpty()) d->m_output_file_parameters << p3;
	if (!p4.isEmpty()) d->m_output_file_parameters << p4;
}

void MSProcessor::setRequiredParameters(const QString &p1, const QString &p2, const QString &p3, const QString &p4)
{
	if (!p1.isEmpty()) d->m_required_parameters << p1;
	if (!p2.isEmpty()) d->m_required_parameters << p2;
	if (!p3.isEmpty()) d->m_required_parameters << p3;
	if (!p4.isEmpty()) d->m_required_parameters << p4;
}

void MSProcessor::setOptionalParameters(const QString &p1, const QString &p2, const QString &p3, const QString &p4)
{
	if (!p1.isEmpty()) d->m_optional_parameters << p1;
	if (!p2.isEmpty()) d->m_optional_parameters << p2;
	if (!p3.isEmpty()) d->m_optional_parameters << p3;
	if (!p4.isEmpty()) d->m_optional_parameters << p4;
}

bool MSProcessor::checkParameters(const QMap<QString, QVariant> &params)
{
	QStringList required;
	required.append(d->m_input_file_parameters);
	required.append(d->m_output_file_parameters);
	required.append(d->m_required_parameters);
	QStringList optional;
	optional.append(d->m_optional_parameters);
	foreach (QString req,required) {
		if (!params.contains(req)) {
			qWarning() << QString("Processor %1 missing required parameter: %2").arg(this->name()).arg(req);
			return false;
		}
	}
	QStringList keys=params.keys();
	foreach (QString key,keys) {
		if ((!required.contains(key))&&(!optional.contains(key))) {
			qWarning() << QString("Processor %1: invalid parameter: %2").arg(this->name()).arg(key);
			return false;
		}
	}
	return true;
}

