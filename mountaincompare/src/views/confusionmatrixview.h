/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/29/2016
*******************************************************/

#ifndef CONFUSIONMATRIXVIEW_H
#define CONFUSIONMATRIXVIEW_H

#include "mvabstractview.h"

#include <mccontext.h>
#include <mvabstractviewfactory.h>

class ConfusionMatrixViewPrivate;
class ConfusionMatrixView : public MVAbstractView {
    Q_OBJECT
public:
    enum PermutationMode {
        NoPermutation,
        RowPermutation,
        ColumnPermutation,
        BothRowBasedPermutation,
        BothColumnBasedPermutation
    };

    friend class ConfusionMatrixViewPrivate;
    ConfusionMatrixView(MVAbstractContext* mvcontext);
    virtual ~ConfusionMatrixView();

    void prepareCalculation() Q_DECL_OVERRIDE;
    void runCalculation() Q_DECL_OVERRIDE;
    void onCalculationFinished() Q_DECL_OVERRIDE;

    MCContext* mcContext();

protected:
    void keyPressEvent(QKeyEvent* evt) Q_DECL_OVERRIDE;
    void prepareMimeData(QMimeData& mimeData, const QPoint& pos) Q_DECL_OVERRIDE;

private slots:
    void slot_permutation_mode_button_clicked();
    void slot_matrix_view_current_element_changed();
    void slot_update_current_elements_based_on_context();
    void slot_export_csv();
    void slot_export_mda();
    void slot_update_selected();

private:
    ConfusionMatrixViewPrivate* d;
};

class ConfusionMatrixViewFactory : public MVAbstractViewFactory {
    Q_OBJECT
public:
    ConfusionMatrixViewFactory(MVMainWindow* mw, QObject* parent = 0);
    QString id() const Q_DECL_OVERRIDE;
    QString name() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    PreferredOpenLocation preferredOpenLocation() const Q_DECL_OVERRIDE;
    MVAbstractView* createView(MVAbstractContext* context) Q_DECL_OVERRIDE;
private slots:
    //void openClipsForTemplate();
};

#endif // CONFUSIONMATRIXVIEW_H
