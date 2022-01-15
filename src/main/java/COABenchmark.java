import org.um.feri.ears.algorithms.Algorithm;
import org.um.feri.ears.algorithms.so.abc.ABC;
import org.um.feri.ears.algorithms.so.cro.CRO;
import org.um.feri.ears.algorithms.so.gwo.GWO;
import org.um.feri.ears.algorithms.so.jade.JADE;
import org.um.feri.ears.algorithms.so.random.RandomWalkAlgorithm;
import org.um.feri.ears.algorithms.so.tlbo.TLBOAlgorithm;
import org.um.feri.ears.benchmark.Benchmark;
import org.um.feri.ears.benchmark.CEC2015Benchmark;
import org.um.feri.ears.benchmark.RPUOed30Benchmark;
import org.um.feri.ears.util.Util;

import java.util.ArrayList;

public class COABenchmark {
    public static void main(String[] args) {
        Util.rnd.setSeed(System.currentTimeMillis()); //set the seed of the random generator
        Benchmark.printInfo = false; //prints one on one results
        //add algorithms to a list
        ArrayList<Algorithm> algorithms = new ArrayList<Algorithm>();
        algorithms.add(new Coyote());
//        algorithms.add(new CRO());
        algorithms.add(new ABC());
//        algorithms.add(new GWO());
//        algorithms.add(new TLBOAlgorithm());
//        algorithms.add(new RandomWalkAlgorithm());
//        algorithms.add(new JADE());

//        RPUOed30Benchmark rpuoed30 = new RPUOed30Benchmark(); // benchmark with prepared tasks and settings
        CEC2015Benchmark cec = new CEC2015Benchmark();

//        rpuoed30.addAlgorithms(algorithms);  // register the algorithms in the benchmark
        cec.addAlgorithms(algorithms);

        cec.run(3);
//        rpuoed30.run(10); //start the tournament with 10 runs/repetitions
    }
}
