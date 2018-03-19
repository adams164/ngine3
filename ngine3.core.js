var Ngine3 = {}
//import Entity from 'ngine3.entity'
Ngine3.Components = {}
Ngine3.Systems = {}
Ngine3.Entity = function Entity () {
    this.id = (+new Date()).toString(16) + (Math.random() * 100000000 | 0).toString(16) + ECS.Entity.prototype._count
    Entity.prototype._count++
    this.components = {}
    return this
}

Ngine3.Entity.prototype._count = 0;

Ngine3.Entity.prototype.addComponent = function addComponent (component){
    this.components[component.name] = component
    return this
}

Ngine3.Entity.prototype.removeComponent = function removeComponent (componentName){
    var name = componentName
    if (typeof componentName === 'function') { 
        name = componentName.prototype.name
    }
    delete this.components[name]
    return this
}
