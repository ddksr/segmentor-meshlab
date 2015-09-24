#ifndef EDIT_SEGMENTOR_H
#define EDIT_SEGMENTOR_H

#include <QObject>
#include <Qt>
#include <QDialog>
#include <QDockWidget>
#include <QUrl>

#include <common/interfaces.h>
#include <meshlab/glarea.h>

#include "segRecoverDialog.h"

class EditSegmentorRecoverPlugin : public QObject, public MeshEditInterface
{
	Q_OBJECT
	Q_INTERFACES(MeshEditInterface)
		
public:
    EditSegmentorRecoverPlugin();
    virtual ~EditSegmentorRecoverPlugin() {}

	static const QString Info();
	
	bool StartEdit(MeshModel &/*m*/, GLArea * /*parent*/);
    void EndEdit(MeshModel &/*m*/, GLArea * /*parent*/);
	void Decorate(MeshModel &, GLArea*);
	
	//void Decorate(MeshModel &/*m*/, GLArea * /*parent*/, QPainter *p) {}
    void mousePressEvent(QMouseEvent *, MeshModel &, GLArea *);
    void mouseMoveEvent(QMouseEvent *, MeshModel &, GLArea *);
    void mouseReleaseEvent(QMouseEvent *, MeshModel &, GLArea *);

	void drawFace(CMeshO::FacePointer fp,MeshModel &m, GLArea *gla, QPainter *p) {}

	// IMPORTANT: because if true, dialog will fail
	bool isSingleMeshEdit() const {return false;}
	
	QFont qFont;

	MeshModel *model; 
	GLArea *gla;
	segRecoverDialog *recoverDialog;
	QSettings *recoverySettings;

	bool editDialogOn;

	QString settingsFile;

};

#endif
