#ifndef EDIT_SEGMENTOR_H
#define EDIT_SEGMENTOR_H

#include <QObject>
#include <common/interfaces.h>
#include <meshlab/glarea.h>
#include <Qt>
#include <QDialog>
#include <QDockWidget>
#include <QUrl>

class EditSegmentorPlugin : public QObject, public MeshEditInterface
{
	Q_OBJECT
	Q_INTERFACES(MeshEditInterface)
		
public:
    EditSegmentorPlugin();
    virtual ~EditSegmentorPlugin() {}

	bool StartEdit(MeshDocument &/*m*/, GLArea * /*parent*/);
	static const QString Info();
	
    void EndEdit(MeshModel &/*m*/, GLArea * /*parent*/) {}
    void Decorate(MeshModel &/*m*/, GLArea * /*parent*/, QPainter *p) {}
    void Decorate (MeshModel &/*m*/, GLArea * ) {}
    void mousePressEvent(QMouseEvent *, MeshModel &, GLArea * ) {}
    void mouseMoveEvent(QMouseEvent *, MeshModel &, GLArea * ) {}
    void mouseReleaseEvent(QMouseEvent *event, MeshModel &/*m*/, GLArea * ) {}
	void drawFace(CMeshO::FacePointer fp,MeshModel &m, GLArea *gla, QPainter *p) {}

	// IMPORTANT: because if true, dialog will fail
	bool isSingleMeshEdit() const {return false;}
	
	QFont qFont;

	MeshDocument *md; 
	GLArea *gla;

private:
	void loadSettings();
	void saveSettings();

	QString settingsFile;

};

#endif
