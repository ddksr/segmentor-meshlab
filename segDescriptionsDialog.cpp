#include <Qt>
#include <QDialog>
#include <QDebug>
#include <vector>
#include <QFileDialog>
#include <QFile>

#include "segDescriptionsDialog.h"
#include "libsegmentor/segmentor.h"

char *modeltypes[] = {"Plane","Superquadric","2nd order surface","Sphere","Cylinder","Cone","Torus","ASQ"};

QString *qmodeltypes[] = {
  new QString("PL"),
  new QString("SQ"),
  new QString("2N"),
  new QString("SP"),
  new QString("CY"),
  new QString("CO"),
  new QString("TO"),
  new QString("ASQ")
};

void segDescriptionsDialog::closeEvent ( QCloseEvent * /*event*/ )
{
  emit closing();
}

segDescriptionsDialog::segDescriptionsDialog(QWidget *parent, GLArea* gl) : QDockWidget(parent) {
  gla = gl;
  segDescriptionsDialog::ui.setupUi(this);

  this->setFeatures(QDockWidget::AllDockWidgetFeatures);
  this->setAllowedAreas(Qt::LeftDockWidgetArea);
  this->setFloating(true);

  QObject::connect(ui.btnGrowStep, SIGNAL(released()), this, SLOT(handleGrowStep()));
  QObject::connect(ui.btnGrowFull, SIGNAL(released()), this, SLOT(handleGrowFull()));
  QObject::connect(ui.btnDelete, SIGNAL(released()), this, SLOT(handleDelete()));
  QObject::connect(ui.btnSelect, SIGNAL(released()), this, SLOT(handleSelect()));

  QObject::connect(ui.btnImport, SIGNAL(released()), this, SLOT(handleImport()));
  QObject::connect(ui.btnExport, SIGNAL(released()), this, SLOT(handleExport()));

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
	for (j = 0; j < seg->descriptions[i].n; j++) {
	  if (k == 0) return;
	  k--;
	}
  }
}

void segDescriptionsDialog::handleSelect() {
  seg->getDrawer()->clear();
  foreach(const QModelIndex &index, 
		  ui.listDescriptions->selectionModel()->selectedIndexes()) {
	int i, j;
	setIJ(index.row(), i, j);
	description *d = seg->descriptions[i].d[seg->descriptions[i].handle[j]];
	seg->getDrawer()->prepare(d->mmodel);
  }
  gla->update();
}

void segDescriptionsDialog::handleDelete() {
  foreach(const QModelIndex &index, 
		  ui.listDescriptions->selectionModel()->selectedIndexes()) {
	int i, j;
	setIJ(index.row(), i, j);
	seg->descriptions[i].delete_description(j);
  }
  fillList();
  reprepareDrawer();
}

void segDescriptionsDialog::handleGrowStep() {
  foreach(const QModelIndex &index, 
		  ui.listDescriptions->selectionModel()->selectedIndexes()) {
	int i, j;
	setIJ(index.row(), i, j);
	description *d = seg->descriptions[i].d[seg->descriptions[i].handle[j]];
	d->grow();
  }
  fillList();
  reprepareDrawer();
}

void segDescriptionsDialog::handleGrowFull() {
  foreach(const QModelIndex &index, 
		  ui.listDescriptions->selectionModel()->selectedIndexes()) {
	int i, j;
	setIJ(index.row(), i, j);
	description *d = seg->descriptions[i].d[seg->descriptions[i].handle[j]];
	while(d->can_grow()) {
	  d->grow();
	}
  }
  fillList();
  reprepareDrawer();
}

void segDescriptionsDialog::handleImport() {
  QString suggestion(".");
  //if(NULL != meshModel)
  //suggestion = PickedPoints::getSuggestedPickedPointsFileName(*meshModel);

  QString filename = QFileDialog::getOpenFileName(this, tr("Load File"), suggestion,
												  "*.segm");
	
  if("" == filename) { return; }
  QFile file(filename);
		
  if(!file.exists()) { return; }

  if(!file.open(QIODevice::ReadOnly)) {
    return;
  }
  
  QTextStream in(&file);
  int paramSize, mtype;
  double value[20];

  while(!in.atEnd()) {
	// close
	QStringList parts = in.readLine().split(" ");
	paramSize = parts.size() - 1;

	mtype = 0;
	while (mtype < NUM_OF_MODELS) {
	  if (qmodeltypes[mtype] == parts.at(0)) { break; }
	  mtype++;
	}
	if (mtype >= NUM_OF_MODELS) {
	  break; // ERROR
	}
	for (int k = 1; k <= paramSize; k++) {
	  value[k] = parts.at(k).toDouble();
	}
	seg->import((MODELTYPE)mtype, value); // TODO: doesn't really work well
  }
  file.close();

  fillList();
}

void segDescriptionsDialog::handleExport() {
  QString suggestion(".");
  // if(NULL != meshModel) {
  // 	suggestion = PickedPoints::getSuggestedPickedPointsFileName(*meshModel); 
  // }

  char *s[20];
  for (int i = 0; i < 20; i++)
	s[i] = new char[30];
  double value[20];

  QString filename = QFileDialog::getSaveFileName(this, tr("Save File"), suggestion, "*.segm");
  qDebug() << filename;
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	return;
  
  QTextStream out(&file);
  
  foreach(const QModelIndex &index, 
		  ui.listDescriptions->selectionModel()->selectedIndexes()) {
	int i, j;
	setIJ(index.row(), i, j);
	description *d = seg->descriptions[i].d[seg->descriptions[i].handle[j]];

	for (int k = 0; k < 20; k++) s[k][0] = 0;
	d->mmodel->parameters(s, value);

	out << *qmodeltypes[(int)d->mmodel->what_model()] << " ";
	for (int k = 0; s[k][0]; k++) {
	  QString param;
	  param.sprintf("%f", value[k]);
	  out << param << " ";
	}
	out << "\n";
  }
}

void segDescriptionsDialog::reprepareDrawer() {
  seg->getDrawer()->clear();
  for (int i = 0; i < seg->numOfDescriptions; i++) {
	for (int j = 0; j < seg->descriptions[i].n; j++) {
	  model* mmodel = seg->descriptions[i].d[seg->descriptions[i].handle[j]]->mmodel;
	  seg->getDrawer()->prepare(mmodel);
	}
  }
  gla->update();
}

// segDescriptionsDialog::~segDescriptionsDialog() {}
