#include "IonChannel.hpp"
#include <gtest/gtest.h>

namespace CSpp {

// Classe dérivée pour les tests nécessitant un accès aux méthodes protégées
class IonChannelTestHelper : public IonChannel {
public:
    using IonChannel::IonChannel;
    using IonChannel::calculateErev;
    using IonChannel::adjustForTemperature;
    using IonChannel::updateGatingVariables;
    using IonChannel::computeCurrent;
    using IonChannel::updateIonConcentrations;
    using IonChannel::getIonConcentrationIn;
    using IonChannel::getIonConcentrationOut;
};

// Test du constructeur et des accesseurs
TEST(IonChannelTest, ConstructorAndAccessors) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    EXPECT_EQ(ionChannel.getGMax(), gMax);
    EXPECT_EQ(ionChannel.getErev(), Erev);
    EXPECT_EQ(ionChannel.getCalciumSensitivity(), calciumSensitivity);
    EXPECT_EQ(ionChannel.getIonConcentrationIn(), ionConcentrationIn);
    EXPECT_EQ(ionChannel.getIonConcentrationOut(), ionConcentrationOut);
}

// Test du calcul du potentiel de réversion (Erev)
TEST(IonChannelTest, CalculateErev) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannelTestHelper ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    double z = 1.0;
    double expectedErev = ionChannel.calculateErev(ionConcentrationIn, ionConcentrationOut, z);

    double R = 8.314;
    double F = 96485.0;
    double T = (temperature + 273.15); // Conversion en Kelvin
    double expectedErevFormula = (R * T) / (z * F) * std::log(ionConcentrationOut / ionConcentrationIn);

    EXPECT_NEAR(expectedErev, expectedErevFormula, 1e-5);
}

// Test de l'ajustement en fonction de la température
TEST(IonChannelTest, AdjustForTemperature) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannelTestHelper ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    double alpha = 1.0;
    double Q10 = 2.0;
    double adjustedAlpha = ionChannel.adjustForTemperature(alpha, 25.0, Q10);

    double expectedAdjustedAlpha = alpha * std::pow(Q10, (25.0 - 20.0) / 10.0);
    EXPECT_NEAR(adjustedAlpha, expectedAdjustedAlpha, 1e-5);
}


// Test du calcul du courant ionique
TEST(IonChannelTest, ComputeCurrent) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    ionChannel.addGatingVariable(0.5, // Valeur initiale
        [](double V) { return 1.0 / (1.0 + std::exp(-(V + 50.0) / 10.0)); },
        [](double V) { return 1.0 / (1.0 + std::exp((V + 50.0) / 10.0)); });

    double V = -60.0;
    double current = ionChannel.computeCurrent(V);

    // Test du courant ionique en fonction de la relation de Nernst
    if (V < Erev) {
        EXPECT_LT(current, 0.0) << "Le courant devrait être négatif pour V = " << V;
    } else {
        EXPECT_GT(current, 0.0) << "Le courant devrait être positif pour V = " << V;
    }
}

// Test de la gestion du gradient ionique
TEST(IonChannelTest, UpdateIonConcentrations) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    double current = 1.0;
    double dt = 0.1;
    ionChannel.updateIonConcentrations(current, dt);

    EXPECT_NE(ionChannel.getIonConcentrationIn(), 140.0);
    EXPECT_NE(ionChannel.getIonConcentrationOut(), 5.0);
}

// Test des cas limites pour les concentrations ioniques
TEST(IonChannelTest, ExtremeIonConcentrations) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 1.0;        // Concentration très faible
    double ionConcentrationOut = 1000.0;    // Concentration très élevée
    ThermalFactors thermalFactors;

    IonChannelTestHelper ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    // Test du calcul du potentiel de réversion dans des conditions extrêmes
    double z = 1.0;
    double extremeErev = ionChannel.calculateErev(ionConcentrationIn, ionConcentrationOut, z);
    EXPECT_NEAR(extremeErev, (8.314 * (temperature + 273.15)) / (z * 96485.0) * std::log(ionConcentrationOut / ionConcentrationIn), 1e-5);
}

// Test des températures extrêmes réalistes
TEST(IonChannelTest, ExtremeTemperatures) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 45.0;  // Température réaliste (45°C)
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannelTestHelper ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    double alpha = 1.0;
    double Q10 = 2.0;
    double adjustedAlpha = ionChannel.adjustForTemperature(alpha, 45.0, Q10);

    double expectedAdjustedAlpha = alpha * std::pow(Q10, (45.0 - 20.0) / 10.0);
    EXPECT_NEAR(adjustedAlpha, expectedAdjustedAlpha, 1e-5);
}

// Test du calcul du courant pour des valeurs extrêmes de potentiel
TEST(IonChannelTest, ComputeCurrentExtremeVoltages) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    ionChannel.addGatingVariable(0.5, // Valeur initiale
        [](double V) { return 1.0 / (1.0 + std::exp(-(V + 50.0) / 10.0)); },
        [](double V) { return 1.0 / (1.0 + std::exp((V + 50.0) / 10.0)); });

    // Test pour un potentiel très positif
    double V_high = 100.0;
    double current_high = ionChannel.computeCurrent(V_high);
    EXPECT_LT(current_high, 0.0) << "Le courant devrait être négatif pour un potentiel très positif";

    // Test pour un potentiel très négatif
    double V_low = -100.0;
    double current_low = ionChannel.computeCurrent(V_low);
    EXPECT_GT(current_low, 0.0) << "Le courant devrait être positif pour un potentiel très négatif";
}

// Test de performance pour le calcul du courant ionique
TEST(IonChannelTest, PerformanceComputeCurrent) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    ionChannel.addGatingVariable(0.5, 
        [](double V) { return 1.0 / (1.0 + std::exp(-(V + 50.0) / 10.0)); },
        [](double V) { return 1.0 / (1.0 + std::exp((V + 50.0) / 10.0)); });

    double V = -60.0;

    // Mesurer le temps avant l'exécution
    auto start = std::chrono::high_resolution_clock::now();

    // Exécution de la méthode computeCurrent à 1 million de reprises pour tester la performance
    for (int i = 0; i < 25000000; ++i) {
        ionChannel.computeCurrent(V);
    }

    // Mesurer le temps après l'exécution
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Afficher le temps d'exécution pour observer les performances
    std::cout << "Temps d'exécution pour 10 million de calculs de courant : " 
              << duration.count() << " secondes." << std::endl;

    // Si le temps dépasse un certain seuil, on échoue le test
    EXPECT_LE(duration.count(), 1.0) << "Le test de performance a échoué. Temps d'exécution trop long.";
}

// Test de performance pour la mise à jour des concentrations ioniques
TEST(IonChannelTest, PerformanceUpdateIonConcentrations) {
    double gMax = 10.0;
    double Erev = -70.0;
    double calciumSensitivity = 0.5;
    double temperature = 37.0;
    double ionConcentrationIn = 140.0;
    double ionConcentrationOut = 5.0;
    ThermalFactors thermalFactors;

    IonChannel ionChannel(gMax, Erev, calciumSensitivity, temperature, ionConcentrationIn, ionConcentrationOut, thermalFactors);

    double current = 1.0;
    double dt = 0.1;

    // Mesurer le temps avant l'exécution
    auto start = std::chrono::high_resolution_clock::now();

    // Exécution de la méthode updateIonConcentrations à 1 million de reprises pour tester la performance
    for (int i = 0; i < 1000000; ++i) {
        ionChannel.updateIonConcentrations(current, dt);
    }

    // Mesurer le temps après l'exécution
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Afficher le temps d'exécution pour observer les performances
    std::cout << "Temps d'exécution pour 1 million de mises à jour de concentrations : " 
              << duration.count() << " secondes." << std::endl;

    // Si le temps dépasse un certain seuil, on échoue le test
    EXPECT_LE(duration.count(), 1.0) << "Le test de performance a échoué. Temps d'exécution trop long.";
}

}
