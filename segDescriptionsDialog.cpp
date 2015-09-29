#include <Qt>
#include <QDialog>
#include <QDebug>
#include <vector>

#include "segDescriptionsDialog.h"
#include "libsegmentor/segmentor.h"

char *modeltypes[] = {"Plane","Superquadric","2nd order surface","Sphere","Cylinder","Cone","Torus"};

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

  stringList = new QStringList();
  listModel = new QStringListModel(*stringList, NULL);
  ui.listDescriptions->setModel(listModel);
  
  fillList();
}

void segDescriptionsDialog::fillList() {
  char *s[20];
  double value[20];
  for (int i = 0; i < 20; i++)
    s[i] = new char[30];
  
  stringList->clear();

  for (int i = 0; i < seg->numOfDescriptions; i++) {
	for (int j = 0; j < seg->descriptions[i].n; j++) {
	  model* mmodel = seg->descriptions[i].d[seg->descriptions[i].handle[j]]->mmodel;
	  QString out = "List " + QString::number(i) + ", desc. " + QString::number(j) + ": ";
	  out = out + modeltypes[(int)mmodel->what_model()] + "(";
	  for (int k = 0; k < 20; k++) s[k][0] = 0;
      seg->descriptions[i].d[seg->descriptions[i].handle[j]]->mmodel->parameters(s,value);
      for (int k = 0; s[k][0]; k++) {
		QString param;
		param.sprintf("%s=%.2f;", s[k], value[k]);
		out += param;
	  }
	  out += ")";
	  stringList->append(out);
	  qDebug() << out;
	}
  }
  
  listModel->setStringList(*stringList);
  for (int i = 0; i < 20; i++)
    delete [] s[i];
}

void segDescriptionsDialog::setIJ(int k, int &i, int &j) {
  for (i = 0; i < seg->numOfDescriptions; i++) {
	for (int j = 0; j < seg->descriptions[i].n; j++) {
	  if (k == 0) return;
	  k--;
	}
  }
}

// segDescriptionsDialog::~segDescriptionsDialog() {}
