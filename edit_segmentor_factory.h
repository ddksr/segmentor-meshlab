#ifndef EditSegmentorFactory_H
#define EditSegmentorFactory_H

#include <QObject>
#include <common/interfaces.h>

class EditSegmentorFactory : public QObject, public MeshEditInterfaceFactory
{
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(MESH_EDIT_INTERFACE_FACTORY_IID)
	Q_INTERFACES(MeshEditInterfaceFactory)

public:
	EditSegmentorFactory();
	virtual ~EditSegmentorFactory() { delete editSegmentor; }

	//gets a list of actions available from this plugin
	virtual QList<QAction *> actions() const;
	
	//get the edit tool for the given action
	virtual MeshEditInterface* getMeshEditInterface(QAction *);
	
	//get the description for the given action
	virtual QString getEditToolDescription(QAction *);
	
private:
	QList <QAction *> actionList;
	QAction *editSegmentor;
};

#endif
