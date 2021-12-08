import org.um.feri.ears.algorithms.Algorithm;
import org.um.feri.ears.algorithms.AlgorithmInfo;
import org.um.feri.ears.algorithms.Author;
import org.um.feri.ears.problems.DoubleSolution;
import org.um.feri.ears.problems.StopCriterionException;
import org.um.feri.ears.problems.Task;
import org.um.feri.ears.util.Comparator.TaskComparator;
import org.um.feri.ears.util.annotation.AlgorithmParameter;
import org.um.feri.ears.util.report.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

import org.um.feri.ears.util.Util;

public class Coyote extends Algorithm {
    @AlgorithmParameter(name = "Number of packs")
    private int numberOfPacks;
    @AlgorithmParameter(name = "Number of Coyotes in a pack")
    private int numberOfCoyotes;
    @AlgorithmParameter(name = "Full population")
    private int popSize;
    @AlgorithmParameter(name = "Probability of Coyote leaving a pack")
    private int leavingProbability;

    Task task;
    private DoubleSolution best; // used to save global best solution
    private ArrayList<DoubleSolution> population;

    public Coyote(int tNumberOfCoyotes, int tNumberOfPacks) {
        numberOfCoyotes = tNumberOfCoyotes;
        numberOfPacks = tNumberOfPacks;
        popSize = numberOfCoyotes * numberOfPacks;
        leavingProbability = (int) (0.005 * (Math.pow(numberOfCoyotes, 2)));

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

        double d = 10;
        double ps = 1/d;
        Vector<Vector<Integer>> lu = new Vector<>();
        Vector<Integer> lower = new Vector<>();
        Vector<Integer> upper = new Vector<>();
        for(int j = 0; j < d; j++){
            lower.add(-10);
            upper.add(10);
        }
        lu.add(lower);
        lu.add(upper);

        var varMin = lower;
        var varMax = upper;
        // Packs init
        var coyotes = new Vector<Vector<Double>>();
        for (int i = 0; i < popSize; i++) {
            var r = new Vector<Double>();
            for (int j = 0; j < varMax.size(); j++) {
                var element = varMin.get(j) + Math.random() * varMax.get(j) - varMin.get(j);
                r.add(element);
            }
            coyotes.add(r);
        }
        var permutation = new Vector<Integer>();
        for (int i = 0; i < popSize; i++) {
            permutation.add(i);
        }
        Collections.shuffle(permutation);
        var ages = new Vector<Integer>();
        var costs = new Vector<Double>();
        var packs = new Vector<Vector<Integer>>();
        var pack = new Vector<Integer>();
        for (Integer integer : permutation) {
            pack.add(integer);
            ages.add(0);
            if (pack.size() % numberOfCoyotes == 0) {
                packs.add(pack);
                pack = new Vector<>();
            }
        }

        // Evaluating coyotes adaptation
        for (int i = 0; i < popSize; i++) {
            costs.add(sphere(coyotes.get(i)));
        }

        var nfEval = popSize;

        // Output vars
        var min = Collections.min(costs);
        var minIndex = costs.indexOf(min);
        var globalParams = coyotes.get(minIndex);

        // Main loop
        var year = 1;
        while (nfEval < 20000) {
            year++;
            // For each pack
            for (int i = 0; i < numberOfPacks; i++) {
                // Coyotes aux
                var coyotesAux = new Vector<Vector<Double>>();
                var newCoyotesAux = new Vector<Vector<Double>>();
                var costsAux = new Vector<Pair<Integer, Double>>();
                var agesAux = new Vector<Integer>();
                var newAgesAux = new Vector<Integer>();
                var tmpPacks = packs.get(i);
                for (int j = 0; j < tmpPacks.size(); j++) {
                    var tmpPack = tmpPacks.get(j);
                    coyotesAux.add(coyotes.get(tmpPack));
                    costsAux.add(new Pair<>(j, costs.get(tmpPack)));
                    agesAux.add(ages.get(tmpPack));
                }
                costsAux.sort(new Comparator<Pair<Integer, Double>>() {
                    @Override
                    public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
                        return o1.getB().compareTo(o2.getB());
                    }
                });
                for (var cost: costsAux) {
                    var index = cost.getA();
                    newCoyotesAux.add(coyotesAux.get(index));
                    newAgesAux.add(agesAux.get(index));
                }
                var cAlpha = newCoyotesAux.get(0);
                // Calculate social tendency
                var tendency = tendency(newCoyotesAux);
                // Update coyotes social condition
                var newCoyotes = new Vector<Vector<Double>>();
                for (int c = 0; c < numberOfCoyotes; c++) {
                    var rc1 = c;
                    while (rc1 == c) {
                        rc1 = Util.nextInt(numberOfCoyotes);
                    }
                    var rc2 = c;
                    while (rc2 == c || rc2 == rc1) {
                        rc2 = Util.nextInt(numberOfCoyotes);
                    }

                    // Try to update the social condition according to alpha and tendency
                    var tmpVector = new Vector<Double>();
                    var firstRandom = Util.nextDouble(0, 1);
                    var secondRandom = Util.nextDouble(0, 1);
                    for (int k = 0; k < newCoyotesAux.get(0).size(); k++) {
                        tmpVector.add(newCoyotesAux.get(c).get(k) + firstRandom * cAlpha.get(k) - newCoyotesAux.get(rc1).get(k) + secondRandom * (tendency.get(k) - newCoyotesAux.get(rc2).get(k)));
                    }
                    // Limit the coyotes
                    tmpVector = limita(tmpVector, (int)d, varMin, varMax);
                    newCoyotes.add(tmpVector);

                    // Evaluate new social condition
                    var newCost = sphere(tmpVector);
                    nfEval++;

                    // Adaptation
                    if (newCost < costsAux.get(c).getB()) {
                        costsAux.setElementAt(new Pair<>(costsAux.get(c).getA(), newCost), c);
                        newCoyotesAux.setElementAt(newCoyotes.get(c), c);
                    }
                }
                // Birth of new coyotes from random parents
                var parentsPermutation = new Vector<Integer>();
                for (int p = 0; p < numberOfCoyotes; p++) {
                    parentsPermutation.add(p);
                }
                Collections.shuffle(parentsPermutation);
                var parents = new Vector<Integer>();
                parents.add(parentsPermutation.get(0));
                parents.add(parentsPermutation.get(1));
                var prob1 = (1-ps)/2;
                var prob2 = prob1;
                var p1 = new Vector<Boolean>();
                var p2 = new Vector<Boolean>();
                var pdr = new Vector<Integer>();
                for (int p = 0; p < d; p++) {
                    pdr.add(p);
                    p1.add(false);
                    p2.add(false);
                }
                Collections.shuffle(pdr);
                p1.setElementAt(true, pdr.get(0));
                p2.setElementAt(true, pdr.get(1));
                var r = new Vector<Double>();
                for (int j = 0; j < d-2; j++) {
                    r.add(Util.nextDouble(0, 1));
                }
                for (int j = 2; j < d; j++) {
                    if (r.get(j-2) < prob1) {
                        p1.setElementAt(true, pdr.get(j));
                    } else {
                        p1.setElementAt(false, pdr.get(j));
                    }
                    if (r.get(j-2) > 1-prob2) {
                        p2.setElementAt(true, pdr.get(j));
                    } else {
                        p2.setElementAt(false, pdr.get(j));
                    }
                }
                // Eventual noise
                var n = new Vector<Boolean>();
                for (int j = 0; j < p1.size(); j++) {
                    var or = p1.get(j) || p2.get(j);
                    n.add(!or);
                }
                // Generate the pup
                var pup = new Vector<Double>();
                var random1 = Util.nextDouble(0, 1);
                for (int j = 0; j < p1.size(); j++) {
                    var element = (p1.get(j) ? 1 : 0) * newCoyotesAux.get(parents.get(0)).get(j) +
                            (p2.get(j) ? 1 : 0) * newCoyotesAux.get(parents.get(1)).get(j) +
                            (n.get(j) ? 1 : 0) * (varMin.get(j) +
                                    Util.nextDouble(0, 1) * (varMax.get(j) - varMin.get(j)));
                    pup.add(element);
                }
                // TODO: Check if the pup will survive
            }
        }

        // the task object holds information about the stopping criterion
        // and information about the problem (number of dimensions, number of constraints, upper and lower bounds...)
        DoubleSolution newSolution;
        best = task.getRandomEvaluatedSolution();  // generate new random solution (number of evaluations is incremented automatically)
        // to evaluate a custom solution or an array use task.eval(mySolution)

        while (!task.isStopCriterion()) { // run until the stopping criterion is not reached
            year++;
            // Execute for each pack
            for (int i = 0; i < numberOfPacks; i++) {

            }
            newSolution = task.getRandomEvaluatedSolution();
            if (task.isFirstBetter(newSolution, best)) { // use method isFirstBetter to check which solution is better (it checks constraints and considers the type of the problem (minimization or maximization))
                best = newSolution;
            }
            task.incrementNumberOfIterations(); // increase the number of generations for each iteration of the main loop
        }
        return best; // return the best solution found
    }

    private double sphere(Vector<Double> coyotes) {
        double sum = 0;
        for (var coyote: coyotes) {
            sum += Math.pow(coyote, 2);
        }
        return sum;
    }

    private Vector<Double> tendency(Vector<Vector<Double>> aux) {
        var medians = new Vector<Double>();
        for (int i = 0; i < aux.get(0).size(); i++) {
            var tmp = new Vector<Double>();
            for(int j = 0; j < aux.size(); j++) {
                tmp.add(aux.get(j).get(i));
            }
            Collections.sort(tmp);
            double median;
            if (tmp.size() % 2 == 0)
                median = (tmp.get(tmp.size() / 2) + tmp.get(tmp.size() / 2 - 1))/2;
            else
                median = tmp.get(tmp.size() / 2);
            medians.add(median);
        }

        return medians;
    }

    private Vector<Double> limita(Vector<Double> tmpVector, int d, Vector<Integer> varMin, Vector<Integer> varMax) {
        var limited = new Vector<Double>();
        for (int i = 0; i < d; i++) {
            limited.add(Math.max(Math.min(tmpVector.get(i), varMax.get(i)), varMin.get(i)));
        }
        return limited;
    }

    private void initPopulation() throws StopCriterionException {
        population = new ArrayList<>();

        for (int i = 0; i < popSize; i++) {
            if (task.isStopCriterion())
                break;
            DoubleSolution newSolution = task.getRandomEvaluatedSolution();
            population.add(newSolution);
        }

        population.sort(new TaskComparator(task));
        best = new DoubleSolution(population.get(0));
    }

    @Override
    public void resetToDefaultsBeforeNewRun() {}
}
