class ConstraintHyperparameters:
    def __init__(self, hyperparameter=1):
        self.hyperparameter = hyperparameter

    def give_hyperparameter(self, hyperparameter):
        self.hyperparameter = hyperparameter


class Hyperparameters:
    def __init__(self, hyperparameters):
        self.hom = ConstraintHyperparameters()
        self.hairpin = ConstraintHyperparameters()
        self.motifGcContent = ConstraintHyperparameters()
        self.keyGcContent = ConstraintHyperparameters()
        self.jointRepeat = ConstraintHyperparameters()
        self.keyInPayload = ConstraintHyperparameters()
        self.similarity = ConstraintHyperparameters()
        self.uniqueJoints = ConstraintHyperparameters()

        for constraint in hyperparameters:
            self.give_hyperparameter(constraint, hyperparameters[constraint])


    def give_hyperparameter(self, constraint, hyperparameter):
        if constraint == 'hom':
                self.hom.give_hyperparameter(hyperparameter)
                return
        if constraint == 'hairpin':
                self.hairpin.give_hyperparameter(hyperparameter)
                return
        if constraint == 'motifGcContent':
                self.motifGcContent.give_hyperparameter(hyperparameter)
                return
        if constraint == 'keyGcContent':
                self.keyGcContent.give_hyperparameter(hyperparameter)
                return
        if constraint == 'similarity':
                self.similarity.give_hyperparameter(hyperparameter)
                return
        if constraint == 'uniqueJoints':
                self.uniqueJoints.give_hyperparameter(hyperparameter)
                return
