Ngine3.Physics = {}
Ngine3.Physics.Components = {}
Ngine3.Physics.Systems = {}

Ngine3.Physics.Components.Position = function Position (pos) {
    pos = pos || [0, 0]
    this.pos = pos
    return this
}
Ngine3.Physics.Components.Position.prototype.name = 'position'

Ngine3.Physics.Components.Motion = function Motion (params) {
    params = params || {}
    this.vel = params.vel || randV()
	this.limit = params.limit || 0
    return this
}
Ngine3.Physics.Components.Motion.prototype.name = 'motion'

Ngine3.Physics.Components.Inertia = function Inertia (params) {
    params = params || {}
    this.force = params.force || randV()
    this.mass = params.mass || 1 + Math.random() * 4
    this.hardness = params.hardness || Math.random()
    return this
}
Ngine3.Physics.Components.Inertia.prototype.name = 'inertia'

Ngine3.Physics.Systems.Physics = {}
Ngine3.Physics.Systems.Physics.run = function physicsRun () {
    var entities = Ngine3.Physics.Systems.Physics.entityList
    var roomWidth = Ngine3.canvas.width
	var roomHeight = Ngine3.canvas.height
	for (var eID in entities) {
        var curEntity = entities[eID]
        var curPos = curEntity.components.position.pos
        var curVel = curEntity.components.motion.vel
        
        var curAccel = timesV(1/ curEntity.components.inertia.mass, curEntity.components.inertia.force)
        
		curVel = sumV(curVel, timesV(this.deltaT, curAccel))
		
        if (curEntity.components.motion.limit && lengthV(curVel) > curEntity.components.motion.limit) {
            curVel = timesV(curEntity.components.motion.limit, normV(curVel))
        }
        
        curPos = sumV(curPos, timesV(this.deltaT, curVel))
		
        
        if (this.boundRoom) {
            if (curPos[0] < curEntity.components.appearance.size) {
                curPos[0] = 2 * curEntity.components.appearance.size - curPos[0]
                curVel[0] = -curVel[0]
            }
            if (curPos[0] > roomWidth - curEntity.components.appearance.size) {
                curPos[0] = 2 * (roomWidth - curEntity.components.appearance.size) - curPos[0]
                curVel[0] = -curVel[0]
            }
            if (curPos[1] < curEntity.components.appearance.size) {
                curPos[1] = 2 * curEntity.components.appearance.size - curPos[1]
                curVel[1] = -curVel[1]
            }
            if (curPos[1] > roomHeight - curEntity.components.appearance.size) {
                curPos[1] = 2 * (roomHeight - curEntity.components.appearance.size) - curPos[1]
                curVel[1] = -curVel[1]
            }
        }
        
        curEntity.components.position.pos = curPos
        curEntity.components.motion.vel = curVel
    }
}

Ngine3.Physics.Systems.Physics.handleCollision = function physicsCollisionHandler (entityPair) {
    var entity1 = Ngine3.Physics.Systems.Physics.entityList[entityPair.id1]
    var entity2 = Ngine3.Physics.Systems.Physics.entityList[entityPair.id2]
    
    if (entity1 && entity2) {
        var vel1 = entity1.components.motion.vel
        var vel2 = entity2.components.motion.vel
        var normal = normV(diffV(entity2.components.position.pos, entity1.components.position.pos))
        var mass1 = entity1.components.inertia.mass
        var mass2 = entity2.components.inertia.mass
        var elasticity = Math.min(1, entity1.components.inertia.hardness * entity2.components.inertia.hardness)

        var normalVel1 = projectV(vel1, normal)
        var normalVel2 = projectV(vel2, normal)
        
        if (dot(diffV(normalVel1,normalVel2),normal) < 0) return
        
        var tangentVel1 = diffV(vel1, normalVel1)
        var tangentVel2 = diffV(vel2, normalVel2)

        var newNormalVel1 = timesV(1/(mass1+mass2), sumV( sumV( timesV(mass1, normalVel1), timesV(mass2, normalVel2)), timesV(mass2 * elasticity, diffV(normalVel2, normalVel1))))
        var newNormalVel2 = timesV(1/(mass1+mass2), sumV( sumV( timesV(mass1, normalVel1), timesV(mass2, normalVel2)), timesV(mass1 * elasticity, diffV(normalVel1, normalVel2))))
        //var newNormalVel1 = timesV(1/(mass1+mass2), sumV( timesV(mass1 - mass2, normalVel1), timesV(2 * mass2, normalVel2)))
        //var newNormalVel2 = timesV(1/(mass1+mass2), sumV( timesV(mass2 - mass1, normalVel2), timesV(2 * mass1, normalVel1)))

        entity1.components.motion.vel = sumV(newNormalVel1, tangentVel1)
        entity2.components.motion.vel = sumV(newNormalVel2, tangentVel2)
    }
}
Ngine3.Physics.Systems.Physics.deltaT = 1
Ngine3.Physics.Systems.Physics.boundRoom = true
Ngine3.Physics.Systems.Physics.entityList = {}