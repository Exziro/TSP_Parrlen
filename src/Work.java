
import sun.nio.ch.ThreadPool;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class Work {
    private int[][] pos=new int[cityNum][2];
    private final static int N = 5; //商旅个数
    private final static int M = 100;// 种群规模
    private final static int min_b = 20;
    private final static int cityNum = 76; // 城市数量，染色体长度 多商旅自己加
    private int T; // 运行代数
    private float[][] distance; // 距离矩阵
    private int[][]distance_break;//断点矩阵
    private float bestDistance; // 最佳长度
    private int[] bestPath; // 最佳路径
    private int[][] oldPopulation; // 父代种群
    private int[][] newPopulation; // 子代种群
    private float[] fitness; // 个体的适应度
    private float[] Pi; // 个体的累积概率
    private float pCorss; // 交叉概率
    private float pMutate; // 变异概率
    private int t;// 当前代数
    private Random random;
    private ExecutorService executorService = Executors.newFixedThreadPool(2);
    private AtomicInteger current=new AtomicInteger(0);
    private AtomicInteger write=new AtomicInteger(0);
    private String name;
    public void readData(String name) {
        File file = new File(name);
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            int line = 1;
            while ((tempString = reader.readLine()) != null) {
                if(line>=7 && !tempString.equals("EOF"))
                {
                    String[] s = tempString.split(" ");
                    // pos1[Integer.parseInt(s[0])-1][0]=Integer.parseInt(s[1]);
                    // pos1[Integer.parseInt(s[0])-1][1]=Integer.parseInt(s[2]);

                    pos[(int)Float.parseFloat(s[0])-1][0]=(int)Float.parseFloat(s[1]);
                    pos[(int)Float.parseFloat(s[0])-1][1]=(int)Float.parseFloat(s[2]);
                }
                line++;
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
        }
    }

    public Work(int t, float corss, float mutate,String name) {
        T = t;
        pCorss = corss;
        pMutate = mutate;
        distance = new float[cityNum][cityNum];
        this.name=name;

        readData(this.name);
        for (int i = 0; i < cityNum - 1; i++) { // 计算距离矩阵, 距离计算方法为伪欧氏距离
            distance[i][i] = 0; // 对角线为0
            for (int j = i + 1; j < cityNum - 1; j++) {
                double rij = Math.sqrt(((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1])));
                distance[i][j] = (int) Math.round(rij); // 四舍五入
                distance[j][i] = distance[i][j];
            }

        }
        //System.out.println("ok");


        distance[cityNum - 1][cityNum - 1] = 0;
        bestDistance = Integer.MIN_VALUE;
        bestPath = new int[cityNum + 1];
        newPopulation = new int[M][cityNum];
        oldPopulation = new int[M][cityNum];
        fitness = new float[M];
        Pi = new float[M];
        random = new Random(System.currentTimeMillis());
    }
    void printgroup(){
        for(int i = 0 ; i < cityNum -1 ; i++){
            for(int j = 0; j < cityNum - 1; j++){
                System.out.print(distance[i][j] + " ");
            }
            System.out.println();
        }
    }
    public void break_ini(){
        //int minb = min_b;
        distance_break = new int[M][N- 1];//多减去一个原因是5个商旅只需要4个断点
        for(int i =0 ; i < M;i++){
            while (true){
                int temp = random.nextInt(min_b - 1);
                distance_break[i][0] = temp;
                for(int j = 1; j <= N-1-1;j++){
                    distance_break[i][j] = random.nextInt(min_b -1) + temp;
                    temp = distance_break[i][j];
                }
                if(cityNum - distance_break[i][N-1-1] <= min_b ){
                    break;
                }
            }

        }
//        for(int i = 0; i < M-1;i++){
//            for(int j = 0; j< N-1;j++){
//                System.out.print(distance_break[i][j]);
//                System.out.print(" ");
//
//            }
//            System.out.println();
//        }


    }


    void initGroup() {
        break_ini();
        int i , j , k;
        ArrayList<Integer> list = new ArrayList<>();
        float dis=0;
        int index=-1;
        for(k = 0; k < M; k++ ){
            int temp = 1;
            for(int a = 0; a < cityNum;a++){
                list.add(a);
            }
            list.remove((Object)oldPopulation[k][0]);
            for(i = 0; i < N ; i++){//断点分段 基因断点虽然是4个 但是实际是五段
                if(i<N-1) {//要确保每次的都是从0开始
                    for (j = temp; j < distance_break[k][i]; j++) {//读取断点后的分段长度
                        dis = Integer.MAX_VALUE;
                        if(j == temp){
                            for (int m = 0; m < list.size(); m++) {
                                if (distance[0][list.get(m)] < dis) {
                                    dis = distance[0][list.get(m)];
                                    index = list.get(m);
                                }
                            }
                        }
                        else{
                            for (int m = 0; m < list.size(); m++) {
                                if (distance[oldPopulation[k][j - 1]][list.get(m)] < dis) {
                                    dis = distance[oldPopulation[k][j - 1]][list.get(m)];
                                    index = list.get(m);
                                }
                            }
                        }
                        oldPopulation[k][j] = index;
                        list.remove((Object) index);
                    }
                    temp = distance_break[k][i];//从下一个起始点开始
                }
                else{
                    for (j = temp; j < (cityNum - temp)+temp; j++) {//读取最后一段分段长度
                        dis = Integer.MAX_VALUE;
                        if(j == temp){
                            for (int m = 0; m < list.size(); m++) {
                                if (distance[0][list.get(m)] < dis) {
                                    dis = distance[0][list.get(m)];
                                    index = list.get(m);
                                }
                            }
                        }
                        else{
                            for (int m = 0; m < list.size(); m++) {
                                if (distance[oldPopulation[k][j - 1]][list.get(m)] < dis) {
                                    dis = distance[oldPopulation[k][j - 1]][list.get(m)];
                                    index = list.get(m);
                                }
                            }
                        }
                        oldPopulation[k][j] = index;
                        list.remove((Object) index);
                    }
                }

            }
        }
        //此处生成的初始内容存在部分重复 需要注意 可采取部分生成 一部分留作他用
        for(int a = 0; a < M ;a++){
            for(int b = 0; b < cityNum;b++){
                System.out.print(oldPopulation[a][b]);
                System.out.print(" ");
            }
            System.out.println();
        }
        /*int i, j, k;
        ArrayList<Integer> list = new ArrayList<>();
        float dis=0;
        int index=-1;

        for (k = 0; k < M; k++) {   // 种群数
            oldPopulation[k][0] = random.nextInt(cityNum);
            if(k<(M)/2)
            {
                //随机生成法
                for (i = 1; i < cityNum; ) {    // 染色体长度
                    oldPopulation[k][i] = random.nextInt(cityNum);
                    for (j = 0; j < i; j++) {
                        if (oldPopulation[k][i] == oldPopulation[k][j]) {
                            break;
                        }
                    }
                    if (j == i) {
                        i++;
                    }
                }
            }else
            {
                //最近领域法
                for(int m=0;m<cityNum;m++)
                {
                    list.add(m);
                }
                list.remove((Object)oldPopulation[k][0]);
                for(i=1;i<cityNum;i++)
                {
                    dis=Integer.MAX_VALUE;
                    for (int m=0;m<list.size();m++) {
                        if(distance[oldPopulation[k][i-1]][list.get(m)]<dis)
                        {
                            dis=distance[oldPopulation[k][i-1]][list.get(m)];
                            index=list.get(m);
                        }
                    }
                    oldPopulation[k][i]=index;
                    list.remove((Object)index);
                }
            }

        }*/
//        int i, j, k;
//        for (k = 0; k < M; k++) {   // 种群数
//            oldPopulation[k][0] = random.nextInt(cityNum);
//            for (i = 1; i < cityNum; ) {    // 染色体长度
//                oldPopulation[k][i] = random.nextInt(cityNum);
//                for (j = 0; j < i; j++) {
//                    if (oldPopulation[k][i] == oldPopulation[k][j]) {
//                        break;
//                    }
//                }
//                if (j == i) {
//                    i++;
//                }
//            }
//        }
    }

    public float evaluate(int[] chromosome,int[] break_point) {//设计计算长度方式
        float len = 0;
        for (int i = 1; i < cityNum; i++) {//无断点长度

            len += distance[chromosome[i - 1]][chromosome[i]];
        }
        len += distance[chromosome[cityNum - 1]][chromosome[0]]; // 回到起点

        for(int i = 0; i < N - 1 ; i++){
            len += distance[chromosome[break_point[i]]][0];
            len += distance[0][chromosome[break_point[i]+1]];
        }//加上断点部分
        return 1/len;
    }

    void countRate() {
        int k;
        double sumFitness = 0;// 适应度总和
        for (k = 0; k < M; k++) {
            sumFitness += fitness[k];
        }
        Pi[0] = (float) (fitness[0] / sumFitness);
        for (k = 1; k < M; k++) {
            Pi[k] = (float) (fitness[k] / sumFitness + Pi[k - 1]);
        }
    }

    public void selectBestChild() {
        int k, i, maxid;
        float maxevaluation;
        maxid = 0;
        maxevaluation = fitness[0];
        for (k = 1; k < M; k++) {
            if (maxevaluation < fitness[k]) {
                maxevaluation = fitness[k];
                maxid = k;
            }
        }
        if (bestDistance < maxevaluation) {
            bestDistance = maxevaluation;
            for (i = 0; i < cityNum; i++) {
                bestPath[i] = oldPopulation[maxid][i];
            }
        }
        System.out.println((int)(1/bestDistance));
        copyGh(0, maxid,newPopulation);// 将当代种群中适应度最高的染色体k复制到新种群中的第一位
    }

    public void selectChild() {
        int k, i, selectId;
        float ran1; // 挑选概率
        for (k = 1; k < M; k++) {
            ran1 = random.nextFloat();
            for (i = 0; i < M; i++) {
                if (ran1 <= Pi[i]) {
                    break;
                }
            }
            selectId = i;
            copyGh(k, selectId,newPopulation);
        }
    }

    public void copyGh(int k, int kk,int[][] target) {
        int i;
        for (i = 0; i < cityNum; i++) {
            target[k][i] = oldPopulation[kk][i];
        }
    }

    public void evolution() {
        selectBestChild();
        selectChild();

        float r;
        for (int k = 1; k < M-1; k+=2) {
            r = random.nextFloat(); // 交叉概率
            if (r < pCorss) { // 交叉
                orderCrossover(k, k+1);
            }
            r = random.nextFloat(); // 变异概率
            if (r < pMutate) {
                variation(k);
            }
            r = random.nextFloat(); // 变异概率
            if (r < pMutate) {
                variation(k+1);
            }
            //余下四6行注释
            SiOrSwap(newPopulation[k],k,0 );
            SiOrSwap(newPopulation[k],k,1 );
            Opt(newPopulation[k],k);
            SiOrSwap(newPopulation[k+1],k+1,1);
            SiOrSwap(newPopulation[k+1],k+1,0 );
            Opt(newPopulation[k+1],k+1);
        }
        search();
        //如果是偶数，会剩下最后一个没有成对，对其经行变异操作
        if(M%2==0)
            variation(M-1);

    }

    //局部搜索算法
    private void search(){
        for(int i=0;i<M;i++)
        {
            executorService.submit(()->{
                Integer old=current.get();
                Integer cold=write.get();
                while (!current.compareAndSet(old,old+1)){
                    old=current.get();
                }
                if(old+1<M)
                {
                    SiOrSwap(newPopulation[old+1],old+1,0 );
                    SiOrSwap(newPopulation[old+1],old+1,1 );
                    Opt(newPopulation[old+1],old+1);
                    Opt(newPopulation[old+1],old+1);
                    while (!write.compareAndSet(cold,cold+1)){
                        cold=write.get();
                    }
                }

                if(old+1==M-1)
                {
                    while (write.get()!=M-1){};
                    synchronized (newPopulation)
                    {
                        newPopulation.notifyAll();
                    }
                }
            });
        }

        while(current.get()!=M)
        {
            synchronized (newPopulation)
            {
                try {
                    newPopulation.wait();
                }catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
        current.set(0);
        write.set(0);
    }

    //opt算法
    public void Opt(int[] chromosome,int index){
        boolean flag=false;
        float fit=fitness[index];
        while (!flag)
        {
            int[] copy=Arrays.copyOfRange(chromosome,0,chromosome.length);
            flag=true;
            for(int i=0;i<cityNum-1;i++)
            {
                for(int j=i+1;j<cityNum;j++)
                {
                    for(int k=0;k<((j-i)%2==0 ? ((j-i)/2):(j-i)/2+1);k++)
                    {
                        int tmp=copy[i+k];
                        copy[i+k]=copy[j-k];
                        copy[j-k]=tmp;
                    }
                    float f=evaluate(copy,distance_break[index]);// 未判定部分，在没确定index是什么之前无法确定
                    if(f>fit)
                    {
                        fit=f;
                        flag=false;
                        chromosome=Arrays.copyOfRange(copy,0,copy.length);

                    }
                }
            }
        }
    }

    //Si或者Swap算法，type==0时为Si算法，type==1时为Swap算法
    public void SiOrSwap(int[] chromosome,int index,int type) {
        ArrayList<Integer> copy  = new ArrayList<>();
        for(int i=0;i<cityNum;i++)
        {
            copy.add(chromosome[i]);
        }
        int index1=0;
        boolean flag=true;
        float fit=fitness[index];
        //表示多少个基因
        for(int i=0;i<cityNum;i++)
        {
            //找出当前基因的位置下标
            for(int j=0;j<cityNum;j++)
            {
                if(chromosome[j]==i)
                {
                    index1=j;
                    break;
                }
            }
            for(int j=0;j<cityNum;j++)
            {
                if(j!=index1)
                {
                    //判断0代表插入，1代表交换
                    if(type==0)
                    {
                        copy.remove((Object)copy.get(index1));
                        copy.add(j,i);
                    }else if(type==1)
                    {
                        copy.set(index1,copy.get(j));
                        copy.set(j,i);
                    }

                    //计算适应度
                    float f=0;
                    for(int m=1;m<cityNum;m++)
                        f+=distance[copy.get(m-1)][copy.get(m)];
                    f+=distance[copy.get(cityNum-1)][0];
                    //调整部分
                    for(int m = 0; m < N -1 ; m++){
                        f += distance[chromosome[distance_break[index][m]]][0];
                        f += distance[0][chromosome[distance_break[index][m]+1]];
                    }
                    f=1/f;

                    //如果适应度大于原来的就替换
                    if(f>fit)
                    {
                        fit=f;
                        for(int m=0;m<cityNum;m++)
                            chromosome[m]=copy.get(m);
                        flag=false;
                    }
                    //还原
                    if(type==0)
                    {
                        copy.remove((Object)copy.get(j));
                        copy.add(index1,i);
                    }
                    else if(type==1)
                    {
                        copy.set(j,copy.get(index1));
                        copy.set(index1,i);
                    }
                }
            }


            if(!flag)
            {
                copy.clear();
                for(int m=0;m<cityNum;m++)
                {
                    copy.add(chromosome[m]);
                }
            }
        }
    }

    public boolean hasElement(int[] a, int b) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] == b) {
                return true;
            }
        }
        return false;
    }


    public void orderCrossover(int k1, int k2) {

        int[] child1 = new int[cityNum];
        int[] child2 = new int[cityNum];
        for(int i=0;i<cityNum;i++)
        {
            child1[i]=-1;
            child2[i]=-1;
        }
        int ran1 = random.nextInt(cityNum);
        int ran2 = random.nextInt(cityNum);
        while (ran1 == ran2) {
            ran2 = random.nextInt(cityNum);
        }
        if (ran1 > ran2) {
            int temp = ran1;
            ran1 = ran2;
            ran2 = temp;
        }
        for (int i = ran1 ;i <= ran2; i++) {  // 生成子代交叉部分
            child1[i]=newPopulation[k1][i];
            child2[i]=newPopulation[k2][i];
        }

        for (int i = 0; i < cityNum; i++) {
            if (i >= ran1 && i <= ran2) {
                continue;
            }
            for (int j = 0; j < cityNum; j++) {
                if (!hasElement(child1,newPopulation[k2][j])) {
                    child1[i]=newPopulation[k2][j];
                    break;
                }
            }

        }

        for (int i = 0; i < cityNum; i++) {
            if (i >= ran1 && i <= ran2) {
                continue;
            }
            for (int j = 0; j < cityNum; j++) {
                if (!hasElement(child2,newPopulation[k1][j])) {
                    child2[i]=newPopulation[k1][j];
                    break;
                }
            }

        }
        newPopulation[k1]=child1;
        newPopulation[k2]=child2;
    }


    /**
     * 随机多次变异
     */
    public void variation(int k) {
        int ran1, ran2, temp;
        int count;
        count = random.nextInt(cityNum); // 变异次数
        for (int i = 0; i < count; i++) {
            ran1 = random.nextInt(cityNum);
            ran2 = random.nextInt(cityNum);
            while (ran1 == ran2) {
                ran2 = random.nextInt(cityNum);
            }
            temp = newPopulation[k][ran1];
            newPopulation[k][ran1] = newPopulation[k][ran2];
            newPopulation[k][ran2] = temp;
        }
    }

    public void run() {
        int i;
        int k;
        // 初始化种群
        initGroup();
        // 计算初始化种群适应度，Fitness[max]
        for (k = 0; k < M; k++) {
            fitness[k] = evaluate(oldPopulation[k],distance_break[k]);
        }
        // 计算初始化种群中各个个体的累积概率，Pi[max]
        countRate();


        System.out.println("init...");
        for (k = 0; k < M; k++) {
            for (i = 0; i < cityNum; i++) {
                System.out.printf("%-3d", oldPopulation[k][i]);
            }
            System.out.println();
        }

        for (t = 0; t < T; t++) {
            evolution();
            // 将新种群newGroup复制到旧种群oldGroup中，准备下一代进化
            //search();
            for (k = 0; k < M; k++) {
                //以下3行注释
                SiOrSwap(newPopulation[k],k,0 );
                SiOrSwap(newPopulation[k],k,1);
                Opt(newPopulation[k],k);
                for (i = 0; i < cityNum; i++) {
                    oldPopulation[k][i] = newPopulation[k][i];
                }
            }
            // 计算种群适应度
            for (k = 0; k < M; k++) {
                fitness[k] = evaluate(oldPopulation[k],distance_break[k]);
            }
            countRate();
        }
        System.out.println("最后种群...");
        for (k = 0; k < M; k++) {
            for (i = 0; i < cityNum; i++) {
                System.out.printf("%-3d", newPopulation[k][i]);
            }
            System.out.println();
        }
        System.out.printf("最小距离：%d\n" , (int)(1/bestDistance));
        System.out.print("最优路径：\n");

        for (i = 0; i < cityNum; i++) {
            System.out.printf("%-3d", bestPath[i]);
            //System.out.println(pos1[bestPath[i]][0]+"\t"+pos1[bestPath[i]][1]);
        }

        executorService.shutdown();
    }
    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        long start = System.currentTimeMillis();

        Work ga = new Work(100, 1, 0.085f,
                "/Users/yixiqi/Downloads/original_data/original_data/pr76.tsp");
        //ga.printgroup();

//        float sum=0;
//        float[] dis=new float[cityNum*cityNum];
//        for(int i=0;i<cityNum;i++)
//        {
//            for(int j=0;j<cityNum;j++)
//            {
//                dis[i*cityNum+j]=ga.distance[i][j];
//            }
//        }
//        Arrays.sort(dis);
//        for(int i=0;i<cityNum;i++)
//        {
//            System.out.println(dis[i+cityNum]);
//            sum+=dis[i+cityNum];
//        }
//        System.out.println(sum);
        ga.run();
//        System.out.println("\n运行:"+(double)(System.currentTimeMillis()-start)/(1000)+"s");

    }

}