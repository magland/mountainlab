#include "sstimeserieswidget.h"
#include <QList>
#include "ssabstractview.h"
#include "sstimeseriesview.h"
#include <QSplitter>
#include <QDebug>
#include <QVBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QProcess>
#include <QInputDialog>
#include <QFileInfo>
#include <QApplication>
#include <QDir>
#include <QTime>
#include <QFileDialog>
#include <QObjectList>
#include <QTextBrowser>

void do_wait(int ms)
{
    QTime timer;
    timer.start();
    while (timer.elapsed() < ms) {
        qApp->processEvents();
    }
}

class SSTimeSeriesWidgetPrivate {
public:
    SSTimeSeriesWidget* q;

    QList<SSAbstractView*> m_views;
    SSAbstractView* m_current_view;
    QSplitter m_splitter;
    QLabel m_info;
    double m_samplerate; //Hz
    QString m_widget_type;
    Mda m_clip_data;

    void set_info(QString txt);
    void connect_view(SSAbstractView* V, bool current_timepoint_only);
    void update_info(SSAbstractView* V);
    void tell_current_view();
};

SSTimeSeriesWidget::SSTimeSeriesWidget(QWidget* parent)
    : QWidget(parent)
{

    d = new SSTimeSeriesWidgetPrivate;
    d->q = this;

    d->m_current_view = 0;
    d->m_samplerate = 0; //eg 20000
    d->m_widget_type = "raw";

    d->m_splitter.setStyleSheet("QSplitter {background: gray;}");
    d->m_splitter.setOrientation(Qt::Vertical);

    d->m_info.setFixedHeight(30);
    d->m_info.setFont(QFont("Arial", 10));
    d->m_info.setMargin(20);

    QWidget* CW = new QWidget();
    QVBoxLayout* CL = new QVBoxLayout();
    CL->setMargin(0);
    CL->setContentsMargins(4, 4, 4, 4);
    CL->setSpacing(0);
    CW->setLayout(CL);

    CL->addWidget(&d->m_splitter);
    CL->addWidget(&d->m_info);

    QVBoxLayout* vlayout = new QVBoxLayout;
    vlayout->setContentsMargins(0, 0, 0, 0);
    vlayout->setSpacing(0);
    setLayout(vlayout);

    setAttribute(Qt::WA_DeleteOnClose);

    vlayout->addWidget(CW);
}

SSTimeSeriesWidget::~SSTimeSeriesWidget()
{
    qDeleteAll(d->m_views);
    delete d;
}

void SSTimeSeriesWidget::addView(SSAbstractView* V)
{

    if (d->m_views.isEmpty())
        V->setTimelineVisible(true);
    else
        V->setTimelineVisible(false);

    d->m_views.append(V);
    d->m_splitter.addWidget(V);
    d->connect_view(V, false);
    V->initialize();

    if (!d->m_current_view) {
        d->m_current_view = V;
    }
    d->tell_current_view();
}

void SSTimeSeriesWidget::setClipData(const Mda& X)
{
    d->m_clip_data = X;
}

SSAbstractView* SSTimeSeriesWidget::view(int index)
{
    if (index >= d->m_views.count())
        return 0;
    if (index < 0)
        return 0;
    return d->m_views[index];
}

void SSTimeSeriesWidgetPrivate::set_info(QString txt)
{
    m_info.setText(txt);
}
void SSTimeSeriesWidgetPrivate::connect_view(SSAbstractView* V, bool current_timepoint_only)
{
    q->connect(V, SIGNAL(currentXChanged()), q, SLOT(slot_current_x_changed()));
    if (!current_timepoint_only) {
        q->connect(V, SIGNAL(currentChannelChanged()), q, SLOT(slot_current_channel_changed()));
        q->connect(V, SIGNAL(xRangeChanged()), q, SLOT(slot_x_range_changed()));
        q->connect(V, SIGNAL(replotNeeded()), q, SLOT(slot_replot_needed()));
        q->connect(V, SIGNAL(selectionRangeChanged()), q, SLOT(slot_selection_range_changed()));
        q->connect(V, SIGNAL(clicked()), q, SLOT(slot_view_clicked()));
    }
}
void SSTimeSeriesWidget::slot_current_x_changed()
{
    SSAbstractView* V = (SSAbstractView*)sender();
    d->update_info(V);

    foreach (SSAbstractView* V2, d->m_views) {
        if (V != V2) {
            if (d->m_views.indexOf(V) >= 0)
                V2->disableSignals();
            V2->setCurrentTimepoint(V->currentTimepoint());
            if (d->m_views.indexOf(V) >= 0)
                V2->enableSignals();
        }
    }
}
void SSTimeSeriesWidget::slot_current_channel_changed()
{
    SSAbstractView* V = (SSAbstractView*)sender();
    d->update_info(V);
}
void SSTimeSeriesWidget::slot_x_range_changed()
{
    SSAbstractView* V = (SSAbstractView*)sender();
    d->update_info(V);

    foreach (SSAbstractView* V2, d->m_views) {
        if (V != V2) {
            V2->disableSignals();
            V2->setXRange(V->xRange());
            V2->enableSignals();
        }
    }
}

void SSTimeSeriesWidget::slot_replot_needed()
{
    /*
	SSAbstractView *V=(SSAbstractView *)sender();

	foreach (SSAbstractView *V2,d->m_views) {
		if (V!=V2) {
			V2->disableSignals();
			V2->setVerticalZoomFactor(V->verticalZoomFactor());
			V2->enableSignals();
		}
	}
	*/
}

void SSTimeSeriesWidget::slot_selection_range_changed()
{
    SSAbstractView* V = (SSAbstractView*)sender();
    d->update_info(V);

    foreach (SSAbstractView* V2, d->m_views) {
        if (V != V2) {
            V2->disableSignals();
            V2->setSelectionRange(V->selectionRange());
            V2->enableSignals();
        }
    }
}

void SSTimeSeriesWidget::slot_view_clicked()
{
    SSAbstractView* W = (SSAbstractView*)sender();
    if (d->m_current_view != W) {
        d->m_current_view = W;
        d->tell_current_view();
    }
}

int get_num_clips_from_process_output(QString output)
{
    int ind1 = output.indexOf("Number of clips: ");
    if (ind1 < 0)
        return 0;
    output = output.mid(ind1 + QString("Number of clips: ").count());
    int ind2 = output.indexOf("\n");
    if (ind2 < 0)
        return 0;
    output = output.mid(0, ind2).trimmed();
    return output.toInt();
}

void SSTimeSeriesWidget::slot_center_on_cursor()
{
    d->m_current_view->centerOnCursor();
}

void SSTimeSeriesWidget::slot_navigation_instructions()
{

    QString html;
    html += "<h3>SPIKESPY</h3>\n";
    html += "<ul>\n";
    html += "<li>Mousewheel scroll: zoom in, zoom out</li>\n";
    html += "<li>+/- keys also zoom in, zoom out</li>\n";
    html += "<li>up/down arrow keys do vertical zoom (on all views together)</li>\n";
    html += "<li>Click and drag to pan</li>\n";
    html += "<li>\"0\" to zoom back to full view</li>\n";
    html += "<li>Left/right arrow keys scroll left/right by one page</li>\n";
    html += "<li>Ctrl + left/right arrow keys move cursor left/right</li>\n";
    html += "<li>Right-click and drag to select a time interval</li>\n";
    html += "<li>ENTER to zoom in to the selected time interval</li>\n";
    //html+="<li>F (when clicked on a view) flips the channel ordering</li>\n";
    html += "</ul>\n";

    QTextBrowser* W = new QTextBrowser;
    W->setAttribute(Qt::WA_DeleteOnClose);
    W->setHtml(html);
    W->show();
    W->resize(800, 800);
    W->move(this->geometry().topLeft() + QPoint(100, 100));
}

void SSTimeSeriesWidgetPrivate::update_info(SSAbstractView* V)
{
    QString str = "";
    Vec2 SR = V->selectionRange();
    if (SR.x >= 0) {
        str = QString("Duration = %1").arg(V->getTimeLabelForX(SR.y - SR.x));
    }
    else {
        double x0 = V->currentX();
        double t0 = V->currentTimepoint();
        if ((x0 >= 0) && (t0 >= 0)) {
            str = QString("Time = %1 (%2)").arg(V->getTimeLabelForX(x0)).arg((long)t0 + 1);

            //str+=QString("; Value = %1  ").arg(V->currentValue());
        }
    }
    double tmp_mb_1 = ((int)jbytesallocated()) * 1.0 / 1000000;
    double tmp_mb_2 = ((int)jnumbytesread()) * 1.0 / 1000000;
    str += QString("       %1/%2 MB used/read").arg((int)tmp_mb_1).arg((int)tmp_mb_2);

    set_info(str);
}

void SSTimeSeriesWidgetPrivate::tell_current_view()
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i] == m_current_view) {
            m_views[i]->setActivated(true);
        }
        else {
            m_views[i]->setActivated(false);
        }
    }
}
