import org.um.feri.ears.algorithms.Algorithm;
import org.um.feri.ears.algorithms.AlgorithmInfo;
import org.um.feri.ears.algorithms.Author;
import org.um.feri.ears.problems.DoubleSolution;
import org.um.feri.ears.problems.StopCriterionException;
import org.um.feri.ears.problems.Task;
import org.um.feri.ears.util.annotation.AlgorithmParameter;
import org.um.feri.ears.util.report.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.um.feri.ears.util.Util;

public class Coyote extends Algorithm {
    @AlgorithmParameter(name = "Number of packs")
    private final int numberOfPacks;
    @AlgorithmParameter(name = "Number of Coyotes in a pack")
    private final int numberOfCoyotes;
    @AlgorithmParameter(name = "Full population")
    private final int popSize;
    @AlgorithmParameter(name = "Probability of Coyote leaving a pack")
    private final double leavingProbability;

    Task task;

    public Coyote(int tNumberOfCoyotes, int tNumberOfPacks) {
        numberOfCoyotes = tNumberOfCoyotes;
        numberOfPacks = tNumberOfPacks;
        popSize = numberOfCoyotes * numberOfPacks;
        leavingProbability = (0.005 * (Math.pow(numberOfCoyotes, 2)));

        ai = new AlgorithmInfo("COA", "Coyote Optimization Algorithm", ""); // add algorithm name
        au = new Author("kenny", "N/A"); // add author information
    }

    public Coyote() {
        this(5, 20);
    }

    @Override
    public DoubleSolution execute(Task task) throws StopCriterionException { // method which executes the algorithm
        this.task = task;
        if (numberOfCoyotes < 3) {
            throw new StopCriterionException("Number of Coyotes per pack should be at least 3");
        }

        double ps = 1 / (double)task.getNumberOfDimensions();
        ArrayList<Double> lower = new ArrayList<>();
        ArrayList<Double> upper = new ArrayList<>();
        for(int j = 0; j < task.getNumberOfDimensions(); j++){
            lower.add(task.getLowerLimit()[j]);
            upper.add(task.getUpperLimit()[j]);
        }

        // Packs init
        ArrayList<DoubleSolution> coyotes = new ArrayList<>();
        for (int i = 0; i < popSize; i++) {
            coyotes.add(task.getRandomEvaluatedSolution());
        }
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < popSize; i++) {
            permutation.add(i);
        }
        Collections.shuffle(permutation);
        ArrayList<Integer> ages = new ArrayList<Integer>();
        ArrayList<Double> costs = new ArrayList<Double>();
        ArrayList<ArrayList<Integer>> packs = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> pack = new ArrayList<Integer>();
        for (Integer integer : permutation) {
            pack.add(integer);
            ages.add(0);
            if (pack.size() % numberOfCoyotes == 0) {
                packs.add(pack);
                pack = new ArrayList<>();
            }
        }

        // Evaluating coyotes adaptation
        for (int i = 0; i < popSize; i++) {
            costs.add(coyotes.get(i).getEval());
        }

        // Output vars
        DoubleSolution newSolution;
        DoubleSolution best = task.getRandomEvaluatedSolution();
        // Main loop
        while (!task.isStopCriterion()) {
            // For each pack
            for (int p = 0; p < numberOfPacks; p++) {
                // Coyotes aux
                ArrayList<DoubleSolution> coyotesAux = new ArrayList<DoubleSolution>();
                ArrayList<DoubleSolution> newCoyotesAux = new ArrayList<DoubleSolution>();
                ArrayList<Pair<Integer, Double>> costsAux = new ArrayList<Pair<Integer, Double>>();
                ArrayList<Integer> agesAux = new ArrayList<Integer>();
                ArrayList<Integer> newAgesAux = new ArrayList<Integer>();
                ArrayList<Integer> tmpPacks = packs.get(p);
                for (int j = 0; j < tmpPacks.size(); j++) {
                    Integer tmpPack = tmpPacks.get(j);
                    coyotesAux.add(coyotes.get(tmpPack));
                    costsAux.add(new Pair<>(j, costs.get(tmpPack)));
                    agesAux.add(ages.get(tmpPack));
                }
                costsAux.sort(Comparator.comparing(Pair::getB));
                for (Pair<Integer, Double> cost : costsAux) {
                    Integer index = cost.getA();
                    newCoyotesAux.add(coyotesAux.get(index));
                    newAgesAux.add(agesAux.get(index));
                }
                // Update coyotes social condition
                ArrayList<DoubleSolution> newCoyotes = new ArrayList<DoubleSolution>();
                DoubleSolution cAlpha = newCoyotesAux.get(0);
                // Calculate social tendency
                ArrayList<Double> tendency = tendency(newCoyotesAux);
                for (int c = 0; c < numberOfCoyotes; c++) {
                    int rc1 = c;
                    while (rc1 == c) {
                        rc1 = Util.nextInt(numberOfCoyotes);
                    }
                    int rc2 = c;
                    while (rc2 == c || rc2 == rc1) {
                        rc2 = Util.nextInt(numberOfCoyotes);
                    }

                    // Try to update the social condition according to alpha and tendency
                    double[] coyoteValues = new double[newCoyotesAux.get(c).getDoubleVariables().length];
                    double firstRandom = Util.nextDouble();
                    double secondRandom = Util.nextDouble();
                    for (int k = 0; k < newCoyotesAux.get(0).getDoubleVariables().length; k++) {
                        coyoteValues[k] = newCoyotesAux.get(c).getDoubleVariables()[k] +
                                firstRandom * cAlpha.getDoubleVariables()[k] - newCoyotesAux.get(rc1).getDoubleVariables()[k] +
                                secondRandom * (tendency.get(k) - newCoyotesAux.get(rc2).getDoubleVariables()[k]);
                    }
                    // Limit the coyotes
                    coyoteValues = task.setFeasible(coyoteValues);
//                    coyoteValues = limita(coyoteValues, task.getNumberOfDimensions(), lower, upper);

                    // Eval new social condition
                    if (task.isStopCriterion()) {
                        break;
                    }
                    newSolution = task.eval(coyoteValues);

                    newCoyotes.add(newSolution);
                    double newCost = newSolution.getEval();

                    // Adaptation
                    if (task.isFirstBetter(newCost, costsAux.get(c).getB())) {
                        costsAux.set(c, new Pair<>(costsAux.get(c).getA(), newCost));
                        newCoyotesAux.set(c, newCoyotes.get(c));
                    }
                }
                // Birth of new coyotes from random parents
                ArrayList<Integer> parents = randomPermutation(numberOfCoyotes);
                double prob = (1 - ps) / 2;
                ArrayList<Boolean> p1 = new ArrayList<Boolean>();
                ArrayList<Boolean> p2 = new ArrayList<Boolean>();
                ArrayList<Integer> pdr = new ArrayList<Integer>();
                for (int i = 0; i < task.getNumberOfDimensions(); i++) {
                    pdr.add(i);
                    p1.add(false);
                    p2.add(false);
                }
                Collections.shuffle(pdr);
                p1.set(pdr.get(0), true);
                p2.set(pdr.get(1), true);
                ArrayList<Double> r = new ArrayList<Double>();
                for (int j = 0; j < task.getNumberOfDimensions() - 2; j++) {
                    r.add(Util.nextDouble(0, 1));
                }
                for (int j = 2; j < task.getNumberOfDimensions(); j++) {
                    if (r.get(j - 2) < prob) {
                        p1.set(pdr.get(j), true);
                    } else {
                        p1.set(pdr.get(j), false);
                    }
                    if (r.get(j - 2) > 1 - prob) {
                        p2.set(pdr.get(j), true);
                    } else {
                        p2.set(pdr.get(j), false);
                    }
                }
                // Eventual noise
                ArrayList<Boolean> n = new ArrayList<Boolean>();
                for (int j = 0; j < p1.size(); j++) {
                    Boolean or = p1.get(j) || p2.get(j);
                    n.add(!or);
                }
                // Generate the pup
                double[] pup = new double[p1.size()];
                for (int j = 0; j < p1.size(); j++) {
                    pup[j] = (p1.get(j) ? 1 : 0) * newCoyotesAux.get(parents.get(0)).getDoubleVariables()[j] +
                            (p2.get(j) ? 1 : 0) * newCoyotesAux.get(parents.get(1)).getDoubleVariables()[j] +
                            (n.get(j) ? 1 : 0) * (lower.get(j) +
                                    Util.nextDouble(0, 1) * (upper.get(j) - lower.get(j)));
                }
                // Check if the pup will survive
                if (task.isStopCriterion()) {
                    break;
                }
                DoubleSolution pupCost = task.eval(pup);

                ArrayList<Integer> worst = new ArrayList<Integer>();
                for (int j = 0; j < costsAux.size(); j++) {
                    if (costsAux.get(j).getB() > pupCost.getEval()) {
                        worst.add(j);
                    }
                }
                if (worst.size() > 0) {
                    ArrayList<Pair<Integer, Integer>> older = new ArrayList<Pair<Integer, Integer>>();
                    for (int j = 0; j < worst.size(); j++) {
                        older.add(new Pair<>(j, newAgesAux.get(j)));
                    }
                    older.sort(Comparator.comparing(Pair::getB));
                    Collections.reverse(older);
                    ArrayList<Integer> which = new ArrayList<Integer>();
                    for (Pair<Integer, Integer> integerIntegerPair : older) {
                        which.add(worst.get(integerIntegerPair.getA()));
                    }
                    newCoyotesAux.set(which.get(0), pupCost);
                    costsAux.set(which.get(0), new Pair<>(which.get(0), pupCost.getEval()));
                    newAgesAux.set(which.get(0), 0);
                }
                // Update the pack info
                for (int j = 0; j < packs.get(p).size(); j++) {
                    Integer packIndex = packs.get(p).get(j);
                    coyotes.set(packIndex, newCoyotesAux.get(j));
                    costs.set(packIndex, costsAux.get(j).getB());
                    ages.set(packIndex, newAgesAux.get(j));
                }
            }
            // Coyote leaving a pack
            if (numberOfPacks > 1) {
                if (Util.nextDouble(0, 1) < leavingProbability) {
                    ArrayList<Integer> rp = randomPermutation(numberOfPacks);
                    Pair<Integer, Integer> rc = new Pair<>(Util.nextInt(numberOfCoyotes), Util.nextInt(numberOfCoyotes));
                    // Coyote leaves the pack and joins a random one
                    Integer aux = packs.get(rp.get(0)).get(rc.getA());
                    Integer tmp = packs.get(rp.get(1)).get(rc.getB());
                    packs.get(rp.get(0)).set(rc.getA(), tmp);
                    packs.get(rp.get(1)).set(rc.getB(), aux);
                }
            }
            // Update ages
            for (int i = 0; i < ages.size(); i++) {
                ages.set(i, ages.get(i) + 1);
            }

            // Output vars
            Double globalMin = Collections.min(costs);
            int ibest = costs.indexOf(globalMin);
            best = coyotes.get(ibest);
            task.incrementNumberOfIterations();
        }
        return best; // return the best solution found
    }

    private ArrayList<Integer> randomPermutation(int num) {
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num; i++) {
            permutation.add(i);
        }
        Collections.shuffle(permutation);
        ArrayList<Integer> parents = new ArrayList<Integer>();
        parents.add(permutation.get(0));
        parents.add(permutation.get(1));
        return parents;
    }

    private ArrayList<Double> tendency(ArrayList<DoubleSolution> pack) {
        ArrayList<Double> medians = new ArrayList<Double>();
        for (int i = 0; i < pack.get(0).getDoubleVariables().length; i++) {
            ArrayList<Double> tmp = new ArrayList<Double>();
            for (DoubleSolution doubleSolution : pack) {
                tmp.add(doubleSolution.getDoubleVariables()[i]);
            }
            Collections.sort(tmp);
            double median;
            if (tmp.size() % 2 == 0)
                median = (tmp.get(tmp.size() / 2) + tmp.get(tmp.size() / 2 - 1)) / 2;
            else
                median = tmp.get(tmp.size() / 2);
            medians.add(median);
        }

        return medians;
    }

    private double[] limita(double[] coyotes, int d, ArrayList<Double> varMin, ArrayList<Double> varMax) {
        double[] limited = new double[coyotes.length];
        for (int i = 0; i < d; i++) {
            limited[i] = Math.max(Math.min(coyotes[i], varMax.get(i)), varMin.get(i));
        }
        return limited;
    }

    @Override
    public void resetToDefaultsBeforeNewRun() {
    }
}
