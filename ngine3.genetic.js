Ngine3.Genetic = {}
Ngine3.Genetic.Components = {}
Ngine3.Genetic.Systems = {}

Ngine3.Genetic.Components.Fitness = function Fitness (fit) {
    this.fit = fit || 1
}
Ngine3.Genetic.Components.Fitness.prototype.name = 'fitness'

Ngine3.Genetic.Systems.GA = {}
Ngine3.Genetic.Systems.GA.selection = function selection (orderedGeneration) {
	var totalFitness = orderedGeneration.reduce((acc, cur) => acc + Ngine3.Genetic.Systems.GA.entityList[cur].components.fitness.fit, 0)
	
	var pickFitness = Math.random() * totalFitness
	var fatherIndex = 0
	for (var fitSoFar = 0; fitSoFar < pickFitness; fatherIndex++) {
		fitSoFar += Ngine3.Genetic.Systems.GA.entityList[orderedGeneration[fatherIndex]].components.fitness.fit
	}
	
	pickFitness = Math.random() * totalFitness
	var motherIndex = 0
	for (var fitSoFar = 0; fitSoFar < pickFitness; motherIndex++) {
		fitSoFar += Ngine3.Genetic.Systems.GA.entityList[orderedGeneration[motherIndex]].components.fitness.fit
	}
	
	return [orderedGeneration[fatherIndex], orderedGeneration[motherIndex]]
}
Ngine3.Genetic.Systems.GA.crossover = function crossover (father, mother) {
	var newGenome = JSON.parse(JSON.stringify(Ngine3.Genetic.Systems.GA.entityList[father].components[genomeComponent]))
	var splicePoint = Math.floor(Math.random() * Object.keys(newGenome).length)
	for (var i = 0; i < splicePoint; i++) {
		currentGene = Object.keys(newGenome)[i]
		newGenome[currentGene] = Ngine3.Genetic.Systems.GA.entityList[mother].components[genomeComponent][currentGene]
	}
	return newGenome
}
Ngine3.Genetic.Systems.GA.mutate = function mutate (genome) {
	for (var key in newGenome){
		if (Math.random() < this.mutationChance) {
			newGenome[key] *= 1 + this.mutationRate * (Math.random() - 0.5)
		}
	}
}
Ngine3.Genetic.Systems.GA.nextGeneration = function nextGeneration () {
	var currentGeneration = Ngine3.Genetic.Systems.GA.entityList
	var curOrderedGeneration = Object.keys(currentGeneration).sort(function (a, b) {
		return currentGeneration[b].components.fitness.fit - currentGeneration[a].components.fitness.fit
	})
	
	var nextGenerationGenomes = []
	
	for (var i = 0; i < this.generationSize; i++) {
		var parents = this.selection(curOrderedGeneration)
		var newGenome = this.crossover(parents[0], parents[1])
		newGenome = this.mutate(newGenome)
		
		nextGenerationGenomes.push(newGenome)
	}
	
	
}
Ngine3.Genetic.Systems.GA.entityList = {}
Ngine3.Genetic.Systems.GA.genomeComponent = ''
Ngine3.Genetic.Systems.GA.generationSize = 0
Ngine3.Genetic.Systems.GA.mutationChance = 0.1
Ngine3.Genetic.Systems.GA.mutationRate = 0.05