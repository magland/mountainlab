/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/17/2016
*******************************************************/

#include "flowlayout.h"
#include "mvcontrolpanel.h"

#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QToolButton>
#include <QVBoxLayout>
#include <QDebug>

class ControlManager {
public:
    void add_group_label(QGridLayout* G, QString label);
    QCheckBox* add_check_box(QGridLayout* G, QString name, QString label, bool val);
    QComboBox* add_combo_box(QGridLayout* G, QString name, QString label);
    QLineEdit* add_int_box(QGridLayout* G, QString name, QString label, int val, int minval, int maxval);
    QLineEdit* add_float_box(QGridLayout* G, QString name, QString label, float val, float minval, float maxval);
    QGroupBox* add_radio_button_group(QGridLayout* G, QString name, QStringList options, QString val);
    QPushButton* add_button(QGridLayout* G, QString name, QString label);
    void add_horizontal_divider_line(QVBoxLayout* layout);

    QVariant get_parameter_value(QString name, const QVariant& defaultval = QVariant());
    void set_parameter_value(QString name, QVariant val);
    void set_parameter_label(QString name, QString text);
    void set_parameter_choices(QString name, QStringList choices);
    void set_parameter_enabled(QString name, bool val);

    QCheckBox* checkbox(QString name);

private:
    QMap<QString, QLineEdit*> m_lineedit_controls;
    QMap<QString, QCheckBox*> m_checkbox_controls;
    QMap<QString, QGroupBox*> m_groupbox_controls;
    QMap<QString, QComboBox*> m_combobox_controls;
    QMap<QString, QPushButton*> m_buttons;
};

class MVControlPanelPrivate {
public:
    MVControlPanel* q;
    ControlManager m_controls;

    QLabel* create_group_label(QString label);
    QAbstractButton* find_action_button(QString name);
};

struct action_button_info {
    QString name;
    QString label;
};
action_button_info abi(QString name, QString label)
{
    action_button_info ret;
    ret.name = name;
    ret.label = label;
    return ret;
}

MVControlPanel::MVControlPanel()
{
    d = new MVControlPanelPrivate;
    d->q = this;

    QFont font = this->font();
    font.setFamily("Arial");
    font.setPixelSize(12);
    this->setFont(font);

    QVBoxLayout* layout = new QVBoxLayout;

    {
        //Views
        layout->addWidget(d->create_group_label("Open Views"));
        FlowLayout* F = new FlowLayout;
        layout->addLayout(F);
        QList<action_button_info> BB;
        BB << abi("open-cluster-details", "Cluster Details");
        BB << abi("open-auto-correlograms", "Auto-Correlograms");
        BB << abi("open-matrix-of-cross-correlograms", "Matrix of Cross-Correlograms");
        BB << abi("open-timeseries-data", "Timeseries Data");
        BB << abi("open-clips", "Clips");
        BB << abi("open-clusters", "Cluster(s)");
        BB << abi("open-firing-events", "Firing Events");
        BB << abi("find-nearby-events", "Find Nearby Events");
        for (int i = 0; i < BB.count(); i++) {
            QToolButton* button = new QToolButton;
            QFont font = button->font();
            font.setPixelSize(14);
            button->setFont(font);
            button->setText(BB[i].label);
            button->setProperty("action_name", BB[i].name);
            connect(button, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            F->addWidget(button);
        }
        d->m_controls.add_horizontal_divider_line(layout);
    }

    {
        //Viewing options
        layout->addWidget(d->create_group_label("Viewing Options"));
        QGridLayout* G = new QGridLayout;
        layout->addLayout(G);

        d->m_controls.add_combo_box(G, "timeseries", "Use timeseries:")->setToolTip("Set the timeseries used for display");
        d->m_controls.add_float_box(G, "cc_max_dt_msec", "Max. dt (ms)", 100, 1, 1e6)->setToolTip("Maximum dt for display of cross-correlograms");
        d->m_controls.add_int_box(G, "clip_size", "Clip size (timepoints)", 150, 1, 1e5)->setToolTip("Set clips size used for display");
        QPushButton* BB = new QPushButton("Update all open views");
        BB->setProperty("action_name", "update_all_open_views");
        QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
        layout->addWidget(BB);

        d->m_controls.add_horizontal_divider_line(layout);
    }

    {
        //Shell splitting
        layout->addWidget(d->create_group_label("Shell Splitting"));
        QGridLayout* G = new QGridLayout;
        layout->addLayout(G);

        d->m_controls.add_check_box(G, "use_shell_split", "Use shell split", false);
        QObject::connect(d->m_controls.checkbox("use_shell_split"), SIGNAL(toggled(bool)), this, SLOT(slot_update_enabled_controls()));
        d->m_controls.add_float_box(G, "shell_increment", "Shell increment", 2, 0.1, 1e6)->setToolTip("Minimum thickness of a peak amplitude shell.");
        d->m_controls.add_int_box(G, "min_per_shell", "Min per shell", 150, 0, 1e6)->setToolTip("Minimum number of points in peak amplitude shell.");
        QPushButton* BB = new QPushButton("Apply Shell Splitting");
        BB->setProperty("action_name", "apply_shell_splitting");
        QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
        layout->addWidget(BB);
        d->m_controls.add_horizontal_divider_line(layout);
    }

    {
        //Event filter
        layout->addWidget(d->create_group_label("Event Filter"));
        QGridLayout* G = new QGridLayout;
        layout->addLayout(G);

        d->m_controls.add_check_box(G, "use_event_filter", "Use event filter", false);
        QObject::connect(d->m_controls.checkbox("use_event_filter"), SIGNAL(toggled(bool)), this, SLOT(slot_update_enabled_controls()));
        d->m_controls.add_float_box(G, "min_detectability_score", "Min detectability score", 0, 0, 1e6)->setToolTip("Filter events by detectability score. Use 0 for no filter.");
        d->m_controls.add_float_box(G, "max_outlier_score", "Max outlier score", 3, 0, 1e6)->setToolTip("Filter events by outlier score. Use 0 for no filter.");
        QPushButton* BB = new QPushButton("Apply Filter");
        BB->setProperty("action_name", "apply_filter");
        QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
        layout->addWidget(BB);
        d->m_controls.add_horizontal_divider_line(layout);
    }

    {
        //Annotate/Merge
        layout->addWidget(d->create_group_label("Annotate / Merge"));

        {
            QPushButton* BB = new QPushButton("Annotate selected (A) ...");
            BB->setProperty("action_name", "annotate_selected");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
        {
            QPushButton* BB = new QPushButton("Merge selected (M)");
            BB->setProperty("action_name", "merge_selected");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
        {
            QPushButton* BB = new QPushButton("Unmerge selected (U)");
            BB->setProperty("action_name", "unmerge_selected");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
    }

    {
        //Export
        layout->addWidget(d->create_group_label("Export"));
        {
            QPushButton* BB = new QPushButton("Export MountainView Document (.mv)...");
            BB->setProperty("action_name", "export_mountainview_document");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
        {
            QPushButton* BB = new QPushButton("Export Original Firings...");
            BB->setProperty("action_name", "export_original_firings");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
        {
            QPushButton* BB = new QPushButton("Export Filtered Firings...");
            BB->setProperty("action_name", "export_filtered_firings");
            QObject::connect(BB, SIGNAL(clicked(bool)), this, SLOT(slot_button_clicked()));
            layout->addWidget(BB);
        }
        d->m_controls.add_horizontal_divider_line(layout);
    }

    slot_update_enabled_controls();

    layout->addStretch(0);
    this->setLayout(layout);
}

MVControlPanel::~MVControlPanel()
{
    delete d;
}

void MVControlPanel::setTimeseriesChoices(const QStringList& names)
{
    d->m_controls.set_parameter_choices("timeseries", names);
}

MVViewOptions MVControlPanel::viewOptions() const
{
    MVViewOptions opts;
    opts.cc_max_dt_msec = d->m_controls.get_parameter_value("cc_max_dt_msec").toDouble();
    opts.clip_size = d->m_controls.get_parameter_value("clip_size").toInt();
    opts.timeseries = d->m_controls.get_parameter_value("timeseries").toString();
    return opts;
}

MVEventFilter MVControlPanel::eventFilter() const
{
    MVEventFilter filter;
    filter.use_event_filter = d->m_controls.get_parameter_value("use_event_filter").toBool();
    filter.max_outlier_score = d->m_controls.get_parameter_value("max_outlier_score").toDouble();
    filter.min_detectability_score = d->m_controls.get_parameter_value("min_detectability_score").toDouble();

    filter.use_shell_split = d->m_controls.get_parameter_value("use_shell_split").toBool();
    filter.shell_increment = d->m_controls.get_parameter_value("shell_increment").toDouble();
    filter.min_per_shell = d->m_controls.get_parameter_value("min_per_shell").toInt();
    return filter;
}

void MVControlPanel::setViewOptions(MVViewOptions opts)
{
    d->m_controls.set_parameter_value("cc_max_dt_msec", opts.cc_max_dt_msec);
    d->m_controls.set_parameter_value("clip_size", opts.clip_size);
    d->m_controls.set_parameter_value("timeseries", opts.timeseries);
}

void MVControlPanel::setEventFilter(MVEventFilter X)
{
    d->m_controls.set_parameter_value("use_event_filter", X.use_event_filter);
    d->m_controls.set_parameter_value("max_outlier_score", X.max_outlier_score);
    d->m_controls.set_parameter_value("min_detectability_score", X.min_detectability_score);
    d->m_controls.set_parameter_value("use_shell_split", X.use_shell_split);
    d->m_controls.set_parameter_value("shell_increment", X.shell_increment);
    d->m_controls.set_parameter_value("min_per_shell", X.min_per_shell);
}

QAbstractButton* MVControlPanel::findButton(const QString& name)
{
    return d->find_action_button(name);
}

void MVControlPanel::slot_update_enabled_controls()
{
    bool use_shell_split = d->m_controls.get_parameter_value("use_shell_split").toBool();
    d->m_controls.set_parameter_enabled("shell_increment", use_shell_split);
    d->m_controls.set_parameter_enabled("min_per_shell", use_shell_split);

    bool use_event_filter = d->m_controls.get_parameter_value("use_event_filter").toBool();
    d->m_controls.set_parameter_enabled("max_outlier_score", use_event_filter);
    d->m_controls.set_parameter_enabled("min_detectability_score", use_event_filter);
}

void MVControlPanel::slot_button_clicked()
{
    QWidget* W = qobject_cast<QWidget*>(sender());
    QString action_name = W->property("action_name").toString();
    if (!action_name.isEmpty()) {
        emit userAction(action_name);
    }
}

void ControlManager::add_group_label(QGridLayout* G, QString label)
{
    int r = G->rowCount();
    QLabel* X = new QLabel(label);
    QFont f = X->font();
    f.setPointSize(18);
    X->setFont(f);
    G->addWidget(X, r, 0, 1, 2);
}

QCheckBox* ControlManager::add_check_box(QGridLayout* G, QString name, QString label, bool val)
{
    int r = G->rowCount();
    QCheckBox* X = new QCheckBox;
    X->setChecked(val);
    X->setText(QString("%1").arg(label));
    G->addWidget(X, r, 1);
    m_checkbox_controls[name] = X;
    X->setProperty("name", name);
    return X;
}

QComboBox* ControlManager::add_combo_box(QGridLayout* G, QString name, QString label)
{
    int r = G->rowCount();
    QComboBox* X = new QComboBox;
    X->setProperty("name", name);
    G->addWidget(new QLabel(label), r, 0);
    G->addWidget(X, r, 1);
    m_combobox_controls[name] = X;
    return X;
}

QLineEdit* ControlManager::add_int_box(QGridLayout* G, QString name, QString label, int val, int minval, int maxval)
{
    Q_UNUSED(minval)
    Q_UNUSED(maxval)
    int r = G->rowCount();
    QLineEdit* X = new QLineEdit;
    X->setText(QString("%1").arg(val));
    G->addWidget(new QLabel(label), r, 0);
    G->addWidget(X, r, 1);
    m_lineedit_controls[name] = X;
    return X;
}

QLineEdit* ControlManager::add_float_box(QGridLayout* G, QString name, QString label, float val, float minval, float maxval)
{
    Q_UNUSED(minval)
    Q_UNUSED(maxval)
    int r = G->rowCount();
    QLineEdit* X = new QLineEdit;
    X->setText(QString("%1").arg(val));
    G->addWidget(new QLabel(label), r, 0);
    G->addWidget(X, r, 1);
    m_lineedit_controls[name] = X;
    return X;
}

QGroupBox* ControlManager::add_radio_button_group(QGridLayout* G, QString name, QStringList options, QString val)
{
    int r = G->rowCount();
    QGroupBox* box = new QGroupBox;
    QHBoxLayout* hlayout = new QHBoxLayout;
    foreach(QString option, options)
    {
        QRadioButton* B = new QRadioButton(option);
        if (option == val)
            B->setChecked(true);
        else
            B->setChecked(false);
        B->setProperty("name", name);
        hlayout->addWidget(B);
    }
    hlayout->addStretch();
    box->setLayout(hlayout);
    m_groupbox_controls[name] = box;
    G->addWidget(box, r, 0, 1, 2);
    return box;
}

QPushButton* ControlManager::add_button(QGridLayout* G, QString name, QString label)
{
    int r = G->rowCount();
    QPushButton* X = new QPushButton(label);
    //X->setFixedHeight(20);
    G->addWidget(X, r, 1);
    X->setProperty("name", name);
    m_buttons[name] = X;
    return X;
}

void ControlManager::add_horizontal_divider_line(QVBoxLayout* layout)
{
    QFrame* line = new QFrame;
    line->setFrameShape(QFrame::HLine); // Horizontal line
    line->setFrameShadow(QFrame::Sunken);
    line->setLineWidth(1);
    //layout->addSpacing(25);
    layout->addWidget(line);
}

QVariant ControlManager::get_parameter_value(QString name, const QVariant& defaultval)
{
    if (m_lineedit_controls.contains(name))
        return m_lineedit_controls[name]->text();
    if (m_checkbox_controls.contains(name))
        return m_checkbox_controls[name]->isChecked();
    if (m_combobox_controls.contains(name))
        return m_combobox_controls[name]->currentText();
    if (m_groupbox_controls.contains(name)) {
        QGroupBox* G = m_groupbox_controls[name];
        QList<QObject*> ch = G->children();
        foreach(QObject * obj, ch)
        {
            QRadioButton* R = dynamic_cast<QRadioButton*>(obj);
            if (R) {
                if (R->isChecked())
                    return R->text();
            }
        }
    }
    return defaultval;
}

void ControlManager::set_parameter_value(QString name, QVariant val)
{
    if (m_lineedit_controls.contains(name))
        return m_lineedit_controls[name]->setText(val.toString());
    if (m_checkbox_controls.contains(name))
        return m_checkbox_controls[name]->setChecked(val.toBool());
    if (m_combobox_controls.contains(name))
        return m_combobox_controls[name]->setCurrentText(val.toString());
    if (m_groupbox_controls.contains(name)) {
        QGroupBox* G = m_groupbox_controls[name];
        QList<QObject*> ch = G->children();
        foreach(QObject * obj, ch)
        {
            QRadioButton* R = dynamic_cast<QRadioButton*>(obj);
            if (R) {
                if (R->text() == val) {
                    R->setChecked(true);
                }
            }
        }
    }
}

void ControlManager::set_parameter_label(QString name, QString text)
{
    if (m_checkbox_controls.contains(name))
        return m_checkbox_controls[name]->setText(text);
    if (m_buttons.contains(name))
        return m_buttons[name]->setText(text);
}

void ControlManager::set_parameter_choices(QString name, QStringList choices)
{
    if (m_combobox_controls.contains(name)) {
        QComboBox* CB = m_combobox_controls[name];
        QString txt = CB->currentText();
        CB->clear();
        foreach(QString choice, choices)
        {
            CB->addItem(choice);
        }
        CB->setCurrentText(txt);
    }
}

void ControlManager::set_parameter_enabled(QString name, bool val)
{
    if (m_lineedit_controls.contains(name))
        return m_lineedit_controls[name]->setEnabled(val);
    if (m_checkbox_controls.contains(name))
        return m_checkbox_controls[name]->setEnabled(val);
    if (m_combobox_controls.contains(name))
        return m_combobox_controls[name]->setEnabled(val);
}

QCheckBox* ControlManager::checkbox(QString name)
{
    return m_checkbox_controls[name];
}

QLabel* MVControlPanelPrivate::create_group_label(QString label)
{
    QLabel* ret = new QLabel(label);
    QFont font = ret->font();
    font.setBold(true);
    font.setFamily("Arial");
    font.setPixelSize(16);
    ret->setFont(font);
    return ret;
}

QAbstractButton* MVControlPanelPrivate::find_action_button(QString name)
{
    QList<QAbstractButton*> buttons = q->findChildren<QAbstractButton*>("", Qt::FindChildrenRecursively);
    foreach(QAbstractButton * B, buttons)
    {
        if (B->property("action_name").toString() == name)
            return B;
    }
    return 0;
}

MVViewOptions MVViewOptions::fromJsonObject(QJsonObject obj)
{
    MVViewOptions ret;
    ret.timeseries = obj["timeseries"].toString();
    if (obj.contains("cc_max_dt_msec"))
        ret.cc_max_dt_msec = obj["cc_max_dt_msec"].toDouble();
    if (obj.contains("clip_size"))
        ret.clip_size = obj["clip_size"].toInt();
    return ret;
}

QJsonObject MVViewOptions::toJsonObject() const
{
    QJsonObject obj;
    obj["timeseries"] = this->timeseries;
    obj["cc_max_dt_msec"] = this->cc_max_dt_msec;
    obj["clip_size"] = this->clip_size;
    return obj;
}

MVEventFilter MVEventFilter::fromJsonObject(QJsonObject obj)
{
    MVEventFilter ret;
    ret.use_shell_split = obj["use_shell_split"].toBool();
    if (obj.contains("shell_increment"))
        ret.shell_increment = obj["shell_increment"].toDouble();
    if (obj.contains("min_per_shell"))
        ret.min_per_shell = obj["min_per_shell"].toInt();
    ret.use_event_filter = obj["use_event_filter"].toBool();
    ret.min_detectability_score = obj["min_detectability_score"].toDouble();
    ret.max_outlier_score = obj["max_outlier_score"].toDouble();
    return ret;
}

QJsonObject MVEventFilter::toJsonObject() const
{
    QJsonObject obj;
    obj["use_shell_split"] = this->use_shell_split;
    obj["shell_increment"] = this->shell_increment;
    obj["min_per_shell"] = this->min_per_shell;
    obj["use_event_filter"] = this->use_event_filter;
    obj["min_detectability_score"] = this->min_detectability_score;
    obj["max_outlier_score"] = this->max_outlier_score;
    return obj;
}
