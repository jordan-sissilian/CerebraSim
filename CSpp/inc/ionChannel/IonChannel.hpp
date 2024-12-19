#pragma once

#include <vector>
#include <map>
#include <functional>
#include <string>
#include <cmath>

namespace CSpp {

// Structure pour gérer les facteurs thermiques (Q10) pour différentes fonctions
struct ThermalFactors {
    double activationQ10;  // Q10 pour les fonctions d'activation
    double inactivationQ10; // Q10 pour les fonctions d'inactivation
    double pumpQ10;        // Q10 pour les pompes ioniques (par exemple Na+/K+)
    
    ThermalFactors(double activation = 2.0, double inactivation = 2.0, double pump = 2.0)
        : activationQ10(activation), inactivationQ10(inactivation), pumpQ10(pump) {}
};

// Structure pour les variables de gating
struct GatingVariable {
    double value;
    std::map<std::string, double> additionalParameters;

    double get(const std::string& propertyName) const {
        auto it = additionalParameters.find(propertyName);
        if (it != additionalParameters.end()) {
            return it->second;
        } else {
            throw std::runtime_error("Property not found: " + propertyName);
        }
    }

    void set(const std::string& propertyName, double value) {
        additionalParameters[propertyName] = value;
    }
};

class IonChannel {
protected:
    double gMax;                // Conductance maximale (mS/cm²)
    double Erev;                // Potentiel de réversion (mV)
    double calciumSensitivity;  // Sensibilité au calcium (facultatif)
    double temperature;         // Température (en °C)
    double ionConcentrationIn;  // Concentration intracellulaire de l'ion (mM)
    double ionConcentrationOut; // Concentration extracellulaire de l'ion (mM)
    double pumpActivity;        // Activité de la pompe ionique (si nécessaire)

    std::vector<GatingVariable> gatingVariables;                        // Variables de gating (m, h, n, etc.)
    std::vector<std::function<double(double)>> alphaFunctions;          // Fonctions alpha pour chaque variable
    std::vector<std::function<double(double)>> betaFunctions;           // Fonctions beta pour chaque variable
    std::vector<std::function<double(double)>> inactivationFunctions;   // Fonctions d'inactivation spécifiques

    ThermalFactors thermalFactors;

    // Fonction pour ajuster les paramètres de gating en fonction de la température et de Q10
    double adjustForTemperature(double alpha, double temperature, double Q10) const;

    // Fonction pour calculer le potentiel de réversion d'un ion (équation de Nernst)
    double calculateErev(double ionConcentrationIn, double ionConcentrationOut, double z) const;

    // Fonction pour gérer l'inactivation (modèle simplifié)
    double inactivation(double voltage, double time, size_t gatingVarIndex) const;

    // Fonction pour gérer l'inactivation à plusieurs états
    double multiStateInactivation(double voltage, size_t gatingVarIndex) const;

public:
    IonChannel(double gMax_, double Erev_, double calciumSensitivity_ = 0.0, double temperature_ = 37.0,
               double ionConcentrationIn_ = 140.0, double ionConcentrationOut_ = 5.0,
               ThermalFactors thermalFactors_ = ThermalFactors());

    // Méthode pour ajouter une variable de gating avec ses fonctions alpha et beta
    void addGatingVariable(double initialValue,
                            std::function<double(double)> alpha,
                            std::function<double(double)> beta,
                            std::function<double(double)> inactivation = nullptr);

    // Mise à jour des variables de gating
    void updateGatingVariables(double V, double dt);

    // Calcul du courant ionique
    double computeCurrent(double V) const;

    // Calcul de la dépendance au calcium
    double calciumDependence(double calciumConcentration) const;

    // Gestion du gradient ionique
    void updateIonConcentrations(double current, double dt);

    double getGMax() const;
    double getErev() const;
    double getCalciumSensitivity() const;
    double getIonConcentrationIn() const;
    double getIonConcentrationOut() const;

    virtual ~IonChannel() = default;
};

} // namespace CSpp
