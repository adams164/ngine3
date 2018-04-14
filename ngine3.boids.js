Ngine3.Boids = {}
Ngine3.Boids.Components = {}
Ngine3.Boids.Systems = {}

Ngine3.Boids.Components.Boid = function Boid (params) {
    params = params || {}
    this.flockWeight = params.flockWeight || .05
    this.alignWeight = params.alignWeight || .125
    this.avoidWeight = params.avoidWeight || 1
    this.goalWeight = params.goalWeight || 0
    this.sightRadius = params.sightRadius * params.sightRadius || 40 * 40
}
Ngine3.Boids.Components.Boid.prototype.name = 'boid'

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
Ngine3.Boids.Systems.Render.entityList = {}

Ngine3.Boids.Systems.Boids = {}
Ngine3.Boids.Systems.Boids.getNeighbors = function getNeighbors (entity) {
    var entities = Ngine3.Boids.Systems.Boids.entityList
    var neighbors = []
    var curPos = entity.components.position.pos
    for (var eID in entities) {
        var distVec = diffV(entities[eID].components.position.pos, curPos)
        if (dot(distVec, distVec) < entity.components.boid.sightRadius && eID != entity.id) {
            neighbors.push(eID)
        }
    }
    return neighbors
}
Ngine3.Boids.Systems.Boids.flockRule = function flockRule (entity, neighbors) {
    var entities = Ngine3.Boids.Systems.Boids.entityList
    var totalPos = neighbors.reduce((accPos, eID) => sumV(accPos, entities[eID].components.position.pos),[0,0])
    var centerPos = timesV(1/neighbors.length,totalPos)
    
    var flockDir = (diffV(centerPos,entity.components.position.pos))
    var flockForce = timesV(entity.components.boid.flockWeight,flockDir)
    
    return flockForce
}
Ngine3.Boids.Systems.Boids.alignRule = function alignRule (entity, neighbors) {
    var entities = Ngine3.Boids.Systems.Boids.entityList
    var totalVel = neighbors.reduce((accVel, eID) => sumV(accVel, entities[eID].components.motion.vel), [0,0])
    var avgVel = timesV(1/neighbors.length,totalVel)
    
    var alignDir = avgVel //diffV(avgVel,entity.components.motion.vel)
    var alignForce = timesV(entity.components.boid.alignWeight,alignDir)
    
    return alignForce
}
Ngine3.Boids.Systems.Boids.avoidRule = function avoidRule (entity, neighbors) {
    var entities = Ngine3.Boids.Systems.Boids.entityList
    var totalPush = [0, 0]
    for (var eID in neighbors) {
        var otherEntity = entities[neighbors[eID]]
        var distVec = diffV(entity.components.position.pos, otherEntity.components.position.pos)
        var dist = lengthV(distVec)
        if (dist < 3 * entity.components.appearance.size) {
            var magAcc = 3 * entity.components.appearance.size - dist
            totalPush = sumV(totalPush, timesV(magAcc, normV(distVec)))
        }
    }
    
    
    var avoidForce = timesV(entity.components.boid.avoidWeight,totalPush)
    
    return avoidForce
}
Ngine3.Boids.Systems.Boids.boundRule = function boundRule (entity) {
    var roomWidth = Ngine3.canvas.width
	var roomHeight = Ngine3.canvas.height
    var pos = entity.components.position.pos
    var retVel = [0,0];
    if (pos[0] < 10){
        retVel[0] = 10;
    }
    else if (pos[0] > (roomWidth - 10)){
        retVel[0] = -10;
    }
    if(pos[1] < 10){
        retVel[1] = 10;
    }
    else if(pos[1] > (roomHeight - 10)){
        retVel[1] = -10;
    }
    return retVel;
}
Ngine3.Boids.Systems.Boids.run = function boidsRun () {
    var entities = Ngine3.Boids.Systems.Boids.entityList
    
    for (var eID in entities) {
        var curEntity = entities[eID]
        var neighbors = this.getNeighbors(curEntity)
        var netForce = this.boundRule(curEntity)
        if (neighbors.length > 0) {
            netForce = sumV(netForce, this.alignRule(curEntity, neighbors))
            netForce = sumV(netForce, this.avoidRule(curEntity, neighbors))
            netForce = sumV(netForce, this.flockRule(curEntity, neighbors))     
        }
        //curEntity.components.inertia.force = netForce
        curEntity.components.motion.vel = sumV(curEntity.components.motion.vel, netForce)
    }
    
}
Ngine3.Boids.Systems.Boids.entityList = {}