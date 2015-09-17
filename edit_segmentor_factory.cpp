#include "edit_segmentor_factory.h"
#include "edit_segmentor_recover.h"

EditSegmentorFactory::EditSegmentorFactory()
{
	editSegmentorRecover = new QAction(QIcon(":/images/segmentor.png"),"Segmentor", this);
	
	actionList << editSegmentorRecover;
	
	foreach(QAction *editAction, actionList)
		editAction->setCheckable(false); 
}

//gets a list of actions available from this plugin
QList<QAction *> EditSegmentorFactory::actions() const
{
	return actionList;
}

//get the edit tool for the given action
MeshEditInterface* EditSegmentorFactory::getMeshEditInterface(QAction *action)
{
	if(action == editSegmentorRecover)
	{
		return new EditSegmentorRecoverPlugin();
	} else assert(0); //should never be asked for an action that isnt here
}

QString EditSegmentorFactory::getEditToolDescription(QAction *)
{
	return EditSegmentorRecoverPlugin::Info();
}

MESHLAB_PLUGIN_NAME_EXPORTER(EditSegmentorFactory)
