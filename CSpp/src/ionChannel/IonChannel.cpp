#include "IonChannel.hpp"

namespace CSpp {

IonChannel::IonChannel(double gMax_, double Erev_, double calciumSensitivity_, double temperature_,
                       double ionConcentrationIn_, double ionConcentrationOut_, ThermalFactors thermalFactors_)
    : gMax(gMax_), Erev(Erev_), calciumSensitivity(calciumSensitivity_), temperature(temperature_),
      ionConcentrationIn(ionConcentrationIn_), ionConcentrationOut(ionConcentrationOut_), pumpActivity(1.0),
      thermalFactors(thermalFactors_) {}

// Fonction pour ajuster les paramètres alpha en fonction de la température et de Q10
double IonChannel::adjustForTemperature(double alpha, double temperature, double Q10) const {
    return alpha * std::pow(Q10, (temperature - 20.0) / 10.0);
}

// Calcul du potentiel de réversion d'un ion selon l'équation de Nernst
double IonChannel::calculateErev(double ionConcentrationIn, double ionConcentrationOut, double z) const {
    double R = 8.314;                   // Constante des gaz parfaits en J/(mol·K)
    double F = 96485.0;                 // Constante de Faraday en C/mol
    double T = (temperature + 273.15);  // Température en Kelvin

    return (R * T) / (z * F) * std::log(ionConcentrationOut / ionConcentrationIn);
}

// Modèle d'inactivation à plusieurs états
double IonChannel::multiStateInactivation(double voltage, size_t gatingVarIndex) const {
    double V_half = -50.0 + gatingVarIndex * 10.0;   // Potentiel de demi-activation dépendant du canal
    double k = 6.0 + gatingVarIndex;                 // Pente de la courbe

    double alpha = 1.0 / (1.0 + std::exp((voltage - V_half) / k)); // Activation
    double beta = 1.0 / (1.0 + std::exp((voltage - V_half) / k));  // Inactivation

    return alpha * (1.0 - beta); // Valeur combinée de l'inactivation
}

// Fonction pour ajouter une nouvelle variable de gating avec ses fonctions alpha et beta
void IonChannel::addGatingVariable(double initialValue, 
                                   std::function<double(double)> alpha, 
                                   std::function<double(double)> beta,
                                   std::function<double(double)> inactivation) {
    GatingVariable newVar;
    newVar.value = initialValue;
    gatingVariables.push_back(newVar);
    alphaFunctions.push_back(alpha);
    betaFunctions.push_back(beta);
    inactivationFunctions.push_back(inactivation);
}

// Mise à jour des variables de gating à chaque pas de temps
void IonChannel::updateGatingVariables(double V, double dt) {
    for (size_t i = 0; i < gatingVariables.size(); ++i) {
        double alpha = alphaFunctions[i](V);
        double beta = betaFunctions[i](V);

        // Ajustement des fonctions alpha et beta en fonction de la température et du Q10
        alpha = adjustForTemperature(alpha, temperature, thermalFactors.activationQ10);
        double tau = 1.0 / (alpha + beta); // Temps de relaxation
        gatingVariables[i].set("tau", tau);
        gatingVariables[i].set("m_inf", alpha / (alpha + beta)); // Valeur d'équilibre

        // Mise à jour de la valeur de la variable de gating
        gatingVariables[i].value += (gatingVariables[i].get("m_inf") - gatingVariables[i].value) * dt / tau;

        // Ajouter l'inactivation pour les canaux qui le nécessitent
        if (inactivationFunctions[i] != nullptr) {
            double inact = inactivationFunctions[i](V);
            gatingVariables[i].value *= inact;  // Moduler la variable de gating par l'inactivation
        } else {
            double multiInact = multiStateInactivation(V, i);
            gatingVariables[i].value *= multiInact;  // Moduler par inactivation multi-états
        }
    }
}

double IonChannel::computeCurrent(double V) const {
    double totalConductance = 0.0;
    for (const auto& gatingVar : gatingVariables) {
        totalConductance += gMax * gatingVar.value; // Somme des conductances
    }

    double ErevDynamic = calculateErev(ionConcentrationIn, ionConcentrationOut, 1.0); // Potentiel de réversion dynamique
    double current = totalConductance * (V - ErevDynamic); 
    if ((V > ErevDynamic && current > 0) || (V < ErevDynamic && current < 0)) {
        current = -current;  // Inverser le signe du courant si la direction est incorrecte
    }

    return current;
}

// Calcul de la dépendance au calcium
double IonChannel::calciumDependence(double calciumConcentration) const {
    return 1.0 / (1.0 + std::exp(-calciumConcentration / 0.1));
}

// Gestion du gradient ionique (mise à jour des concentrations ioniques)
void IonChannel::updateIonConcentrations(double current, double dt) {
    // Modèle simplifié de pompe ionique Na+/K+-ATPase pour réguler les concentrations
    double pumpEfficiency = 0.5 * std::pow(thermalFactors.pumpQ10, (temperature - 20.0) / 10.0); // Température ajustée pour la pompe
    ionConcentrationIn += (current * dt * pumpEfficiency);                                       // Modification de la concentration intracellulaire
    ionConcentrationOut -= (current * dt * pumpEfficiency);                                      // Modification de la concentration extracellulaire
}

double IonChannel::getGMax() const { return gMax; }
double IonChannel::getErev() const { return Erev; }
double IonChannel::getCalciumSensitivity() const { return calciumSensitivity; }
double IonChannel::getIonConcentrationIn() const { return ionConcentrationIn; }
double IonChannel::getIonConcentrationOut() const { return ionConcentrationOut; }

}