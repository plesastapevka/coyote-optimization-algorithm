import org.um.feri.ears.algorithms.Algorithm;
import org.um.feri.ears.algorithms.AlgorithmInfo;
import org.um.feri.ears.algorithms.Author;
import org.um.feri.ears.problems.DoubleSolution;
import org.um.feri.ears.problems.StopCriterionException;
import org.um.feri.ears.problems.Task;
import org.um.feri.ears.util.Comparator.TaskComparator;
import org.um.feri.ears.util.annotation.AlgorithmParameter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Vector;

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

        int d = 10;
        Vector<Vector<Integer>> lu = new Vector<Vector<Integer>>();
        Vector<Integer> lower=new Vector<Integer>();
        Vector<Integer> upper=new Vector<Integer>();
        for(int j = 0; j < d; j++){
            lower.add(-10);
            upper.add(10);
        }
        lu.add(lower);
        lu.add(upper);

        var varMin = lower;
        var varMax = upper;

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
        var packs = new Vector<Vector<Integer>>();
        var pack = new Vector<Integer>();
        for (int i = 0; i < permutation.size(); i++) {
            if (i % numberOfCoyotes == 0) {
                packs.add(pack);
                pack = new Vector<Integer>();
            }
            pack.add(permutation.get(i));
        }

        // the task object holds information about the stopping criterion
        // and information about the problem (number of dimensions, number of constraints, upper and lower bounds...)
        DoubleSolution newSolution;
        best = task.getRandomEvaluatedSolution();  // generate new random solution (number of evaluations is incremented automatically)
        // to evaluate a custom solution or an array use task.eval(mySolution)
        int year = 1;
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
