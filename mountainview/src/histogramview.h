/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef HISTOGRAMVIEW_H
#define HISTOGRAMVIEW_H

#include <QWidget>

/** \class HistogramView
 *  \brief View of a histogram from an arbitrary data array
 */
class HistogramViewPrivate;
class HistogramView : public QWidget {
    Q_OBJECT
public:
    friend class HistogramViewPrivate;
    explicit HistogramView(QWidget* parent = 0);
    virtual ~HistogramView();

    ///The data to view
    void setData(const QList<float> values);
    ///Alternative specification
    void setData(int N, float* values);
    ///Set evenly spaced bins
    void setBins(float bin_min, float bin_max, int num_bins);
    ///Auto set evenly spaced bins based on range of data (call setData first)
    void autoSetBins(int num_bins);
    ///The color for filling the histogram bars
    void setFillColor(const QColor& col);
    ///The edge color for the bars
    void setLineColor(const QColor& col);
    ///Displayed title of histogram
    void setTitle(const QString& title);
    ///Controls background and highlighting colors. For consistent app look. TODO: Low-priority handle this in appropriate way (setColors)
    void setColors(const QMap<QString, QColor>& colors);

    ///Set this as the current histogram (affects highlighting)
    void setCurrent(bool val);
    ///Set this as among the selected histograms (affects highlighting)
    void setSelected(bool val);

    QImage renderImage(int W, int H);

protected:
    void paintEvent(QPaintEvent* evt);
    void mousePressEvent(QMouseEvent* evt);
    void mouseMoveEvent(QMouseEvent* evt);
    void enterEvent(QEvent* evt);
    void leaveEvent(QEvent* evt);
    void mouseDoubleClickEvent(QMouseEvent* evt);
signals:
    ///The histogram was somewhere clicked
    void clicked();
    ///The histogram was somewhere control-clicked
    void control_clicked();
    ///The histogram was somewhere double-clicked (or enter was pressed, not sure??)
    void activated();

    void signalExportHistogramMatrixImage();

private slots:
    void slot_context_menu(const QPoint& pos);

private:
    HistogramViewPrivate* d;
};

#endif // HISTOGRAMVIEW_H
