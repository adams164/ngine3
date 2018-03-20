Ngine3.Boids = {}
Ngine3.Boids.Components = {}
Ngine3.Boids.Systems = {}

Ngine3.Boids.Components.Position = function Position (pos) {
    pos = pos || [0, 0]
    this.pos = pos
    return this
}
Ngine3.Boids.Components.Position.prototype.name = 'position'

Ngine3.Boids.Components.Physics = function Physics (params) {
    params = params || {}
    this.vel = params.vel || randV()
    this.mass = params.mass || 1 + Math.random() * 4
    this.hardness = params.hardness || Math.random()
    return this
}
Ngine3.Boids.Components.Physics.prototype.name = 'physics'

Ngine3.Boids.Components.Boid = function Boid (params) {
    params = params || {}
    this.flockWeight = params.flockWeight || 5
    this.alignWeight = params.alignWeight || 1
    this.avoidWeight = params.avoidWeight || 1
}

Ngine3.Boids.Components.Appearance = function Appearance (params) {
    params = params || {}
    this.color = params.color || {r: 0, g: 0, b: 0}
    this.size = params.size || (Math.random() * 5) + 5
    return this
}
Ngine3.Boids.Components.Appearance.prototype.name = 'appearance'

Ngine3.Boids.Systems.Render = {}
Ngine3.Boids.Systems.Render.run = function renderRun () {
    var entities = Ngine3.Boids.Systems.Render.entityList
    
    Ngine3.context.clearRect(0, 0, Ngine3.canvas.width, Ngine3.canvas.height)
    
    for (var eID in entities) {
        var curEntity = entities[eID]
        var fillStyle = 'rgb(' + [
            curEntity.components.appearance.color.r,
            curEntity.components.appearance.color.g,
            curEntity.components.appearance.color.b
        ] + ')'
        
        Ngine3.context.fillStyle = fillStyle
        
        Ngine3.context.beginPath()
        Ngine3.context.arc(curEntity.components.position.pos[0], curEntity.components.position.pos[1], curEntity.components.appearance.size, 0, 2 * Math.PI)
        Ngine3.context.fill()
    }
}
Ngine3.Boids.Systems.Render.entityList={}

Ngine3.Boids.Systems.Physics = {}
Ngine3.Boids.Systems.Physics.run = function physicsRun () {
    var entities = Ngine3.Boids.Systems.Physics.entityList
    var roomWidth = Ngine3.canvas.width;
	var roomHeight = Ngine3.canvas.height;
	for (var eID in entities) {
        var curEntity = entities[eID]
        var curPos = curEntity.components.position.pos
        var curVel = curEntity.components.physics.vel
        curPos = sumV(curPos, timesV(Ngine3.Boids.Systems.Physics.deltaT, curVel))
		
		
        curVel = sumV(curVel, timesV(Ngine3.Boids.Systems.Physics.deltaT, [0, 0.01]))
        
		
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
        
        curEntity.components.position.pos = curPos
        curEntity.components.physics.vel = curVel
    }
}

Ngine3.Boids.Systems.Physics.handleCollision = function physicsCollisionHandler (entityPair) {
    var entity1 = Ngine3.Boids.Systems.Physics.entityList[entityPair.id1]
    var entity2 = Ngine3.Boids.Systems.Physics.entityList[entityPair.id2]
    
    if (entity1 && entity2) {
        var vel1 = entity1.components.physics.vel
        var vel2 = entity2.components.physics.vel
        var normal = normV(diffV(entity2.components.position.pos, entity1.components.position.pos))
        var mass1 = entity1.components.physics.mass
        var mass2 = entity2.components.physics.mass
        var elasticity = Math.min(1, entity1.components.physics.hardness * entity2.components.physics.hardness)

        var normalVel1 = projectV(vel1, normal)
        var normalVel2 = projectV(vel2, normal)
        
        if (dot(diffV(normalVel1,normalVel2),normal) < 0) return
        
        var tangentVel1 = diffV(vel1, normalVel1)
        var tangentVel2 = diffV(vel2, normalVel2)

        var newNormalVel1 = timesV(1/(mass1+mass2), sumV( sumV( timesV(mass1, normalVel1), timesV(mass2, normalVel2)), timesV(mass2 * elasticity, diffV(normalVel2, normalVel1))))
        var newNormalVel2 = timesV(1/(mass1+mass2), sumV( sumV( timesV(mass1, normalVel1), timesV(mass2, normalVel2)), timesV(mass1 * elasticity, diffV(normalVel1, normalVel2))))
        //var newNormalVel1 = timesV(1/(mass1+mass2), sumV( timesV(mass1 - mass2, normalVel1), timesV(2 * mass2, normalVel2)))
        //var newNormalVel2 = timesV(1/(mass1+mass2), sumV( timesV(mass2 - mass1, normalVel2), timesV(2 * mass1, normalVel1)))

        entity1.components.physics.vel = sumV(newNormalVel1, tangentVel1)
        entity2.components.physics.vel = sumV(newNormalVel2, tangentVel2)
    }
}
Ngine3.Boids.Systems.Physics.deltaT = 1
Ngine3.Boids.Systems.Physics.entityList = {}


Ngine3.Boids.Systems.Collision = {}
Ngine3.Boids.Systems.Collision.run = function collisionRun () {
    var entities = Ngine3.Boids.Systems.Collision.entityList
    var checkedIDs = {}
    var collisionPairs = []
    for (var eID in entities) {
        var curEntity = entities[eID]
        checkedIDs[eID] = true
        for (var eID2 in entities) {
            if (!checkedIDs[eID2]) {
                var entity2 = entities[eID2]
                var distVec = diffV(curEntity.components.position.pos, entity2.components.position.pos)
                var sizes = curEntity.components.appearance.size + entity2.components.appearance.size
                if (dot(distVec, distVec) < sizes * sizes) {
                    collisionPairs.push({id1: eID, id2: eID2})
                }
            }
        }
    }
    
    for (var i in collisionPairs) {
        Ngine3.Boids.Systems.Collision.triggerEvent(collisionPairs[i])
    }
}
Ngine3.Boids.Systems.Collision.triggerEvent = function collisionEvent (entityPair) {
    for (var i in Ngine3.Boids.Systems.Collision.handlers) {
        Ngine3.Boids.Systems.Collision.handlers[i](entityPair)
    }
}
Ngine3.Boids.Systems.Collision.handlers = []
Ngine3.Boids.Systems.Collision.entityList = {}