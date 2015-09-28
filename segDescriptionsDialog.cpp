#include <Qt>
#include <QDialog>

#include "segDescriptionsDialog.h"
#include "libsegmentor/segmentor.h"

void segDescriptionsDialog::closeEvent ( QCloseEvent * /*event*/ )
{
  emit closing();
}

segDescriptionsDialog::segDescriptionsDialog(QWidget *parent) : QDockWidget(parent) {
  segDescriptionsDialog::ui.setupUi(this);

  this->setFeatures(QDockWidget::AllDockWidgetFeatures);
  this->setAllowedAreas(Qt::LeftDockWidgetArea);
  this->setFloating(true);

  seg = Segmentor::Instance();

  
}

void segDescriptionsDialog::fillList() {
  
}

// segDescriptionsDialog::~segDescriptionsDialog() {}
