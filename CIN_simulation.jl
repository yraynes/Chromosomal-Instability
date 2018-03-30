#!/contrib/projects/julia/julia
using Distributions

type Lineage
  fitness::Float64
  size::Int64
  state::Vector{Float64}
  msi_rate::Float64
  cin_rate::Float64
  ben_count::Int64
end

function average_fitness(population::Dict{Array{Float64,1}, Lineage})
  popN = 0
  popW = 0.0
  varianceW = 0.0
  for k in values(population)
    popN+=k.size
    popW+=k.size * k.fitness
  end
  popW = popW/popN
  for k in values(population)
    varianceW+=(k.fitness-popW)^2*k.size/popN
  end
  return popW, popN, varianceW
end

function cancer_state(population::Dict{Array{Float64,1}, Lineage}, threshold::Int64)
  popN, cancerN = 0,0
  # winner = Float64
  # win_freq = Float64
  # carcinoma = Dict{Array{Float64,1}, Lineage}()
  for (key,line) in population
    popN+=line.size
    if line.ben_count >= threshold
      cancerN+=line.size
      # carcinoma[key] = Lineage(line.fitness, line.size, line.state, line.msi_rate, line.cin_rate, line.ben_count)
    end
  end
  # winner, win_freq = find_winner(carcinoma, cancerN)
  cancerF = Float64
  cancerF = cancerN/popN
  return cancerF#, winner, win_freq
end

function find_winner(population::Dict{Array{Float64,1}, Lineage}, N::Int64)
  mutations = Dict{Float64, Int64}()
  for line in values(population)
    if line.state[3] in keys(mutations)
      mutations[line.state[3]]+=line.size
    else
      mutations[line.state[3]] = line.size
    end
  end

  winner = Int64
  size = 0
  for (key, value) in mutations
    if value > size
      size = value
      winner = key
    end
  end
  return winner, size/N
end

function assay_mutation_rate(population::Dict{Array{Float64,1}, Lineage})
  popN, msiN, cinN = 0,0,0
  for line in values(population)
    popN+=line.size
    if line.msi_rate>1 msiN+=line.size end
    if line.cin_rate>0 cinN+=line.size end
  end
  msiF = msiN/popN
  cinF = cinN/popN
  return msiF, cinF
end

#@iprofile begin
function add_to_pop(population::Dict{Array{Float64,1}, Lineage}, newstate::Vector{Float64}, newfit::Float64, count::Int64, ben_count::Int64)
  if newstate in keys(population)
    population[newstate].size += count
  else
    population[newstate] = Lineage(newfit, count, newstate, newstate[4], newstate[5], ben_count)
  end
  return population
end
#end

function make_segment(ben_dis::Int64, cin_mutation_count::Int64, del_dis::Int64)
  segments = zeros(Int64, cin_mutation_count,2)
#  ben_count = Int64[]
  del_count = Int64[]

  if ben_dis == 1
    ben_count = fill(1, cin_mutation_count)
  elseif ben_dis == 2
    ben_count = rand(Poisson(1),cin_mutation_count)
  elseif ben_dis == 3
    ben_count = rand(Geometric(1/2),cin_mutation_count)
  elseif ben_dis == 4
    ben_count = rand(BetaBinomial(50, .1, 4.9),cin_mutation_count)
  end

  if del_dis == 1
      del_count = fill(100, cin_mutation_count)
  elseif del_dis == 2
      del_count = rand(Poisson(100),cin_mutation_count)
  elseif del_dis == 3
      del_count = rand(Geometric(1/100),cin_mutation_count)
  elseif del_dis == 4
      del_count = rand(BetaBinomial(500, .1, .4),cin_mutation_count)
  end

  for i = 1:cin_mutation_count
    segments[i] = ben_count[i]
    segments[i+cin_mutation_count] = del_count[i]
  end
  return segments
end

# @iprofile begin
function mutate_population(population::Dict{Array{Float64,1}, Lineage}, Utot::Array{Float64,1}, ben_dis::Int64, del_dis::Int64, sb::Float64, sd::Float64, msi_multiplier::Float64, cin_multiplier::Float64)
  new_population = Dict{Array{Float64,1}, Lineage}()
  for(key,line) in population
    mutations = rand(Multinomial(line.size, [line.msi_rate*Utot[1],line.msi_rate*Utot[2],line.msi_rate*Utot[3], line.msi_rate*Utot[4], line.msi_rate*Utot[5], 1-line.msi_rate*sum(Utot)]))
    #mutations = [driver1 count, deleterious count, driver2 count, msi count, cin count, non-mutated count]

    if mutations[1] > 0
      #assign driver1 mutations

      newstate = copy(line.state)
      newstate[1] += sb
      new_population = add_to_pop(new_population, newstate, line.fitness + sb, mutations[1], line.ben_count+1)
    end

      #assign deleterious mutations
    if mutations[2] > 0
      newstate = copy(line.state)
      newstate[2] += -sd
      new_population = add_to_pop(new_population, newstate, line.fitness - sd, mutations[2], line.ben_count)
    end

      #assign driver2 mutations
      #could be used for a different kind of driver, i.e. oncogenes vs tum. sups.
#     if mutations[3] > 0
#       newstate = copy(line.state)#
#       newstate[3] += sb
#       newfit = line.fitness + sb
#       new_population = add_to_pop(new_population, newstate, newfit, mutations[3], line.ben_count+1)
#     end

    if mutations[4] > 0
        #assign MSI mutations
      newstate = copy(line.state)#
      newstate[4] = msi_multiplier
      new_population = add_to_pop(new_population, newstate, line.fitness, mutations[4], line.ben_count)
    end

      #assign CIN mutations
    if mutations[5] > 0
      newstate = copy(line.state)#
      newstate[5] = cin_multiplier
      new_population = add_to_pop(new_population, newstate, line.fitness, mutations[5], line.ben_count)
    end

    if line.cin_rate > 0
      cin_mutation_count = rand(Poisson(line.cin_rate * line.size))
      segments = make_segment(ben_dis, cin_mutation_count, del_dis)
      for s = 1:cin_mutation_count
        newstate = copy(line.state)
        newstate[1]+=segments[s,:][1] * sb
        newstate[2]-=segments[s,:][2] * sd
        newval = segments[s,:][1] * sb- segments[s,:][2] * sd
        if newval > newstate[3]
          newstate[3] = newval
        end
        newfit = line.fitness +segments[s,:][1] * sb- segments[s,:][2] * sd
        new_population = add_to_pop(new_population, newstate, newfit, 1, line.ben_count+segments[s,:][1])
      end
    else
      cin_mutation_count = 0
    end

    new_size = mutations[6] - cin_mutation_count
    if  new_size > 0
      new_population = add_to_pop(new_population, line.state, line.fitness, new_size, line.ben_count)
    end
  end
  return new_population
end
# end

function multinomial_rand_samp(n::Int, p::Vector{Float64})
    k = length(p)
    x = Int64[]
    rp = 1.0  # remaining total probability
    i = 0
    km1 = k - 1
    tmp = Float64[]
    while i < km1 && n > 0
        i += 1
        pi = p[i]
        if pi < rp
            xi = rand(Binomial(n, max(pi / rp,0)))
            push!(x, xi)
            push!(tmp, pi)
            n -= xi
            rp -= pi
        else
            push!(x, n)
            n = 0
            rp = 0.0
        end
    end

    if i == km1
        push!(x,n)
    else
        for j = i+1 : k
        push!(x,0)
        end
    end

    return x
end

#@iprofile begin
function wright_fisher_reproduction(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[]
  #key_list = String[]

  for(key,line) in population
    push!(proby_list, line.size/N0 * line.fitness/popw)
    push!(lineage_list, line)
  end
  new_counts_list = multinomial_rand_samp(Nnew, proby_list)
  new_population = Dict{Array{Float64,1}, Lineage}()

  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].msi_rate, lineage_list[i].cin_rate, lineage_list[i].ben_count)
    end
  end
  return new_population
end
#end

#@iprofile begin
function simulate(ben_dis::Int64, del_dis::Int64)
  #
  const popN0 = 1000000
  const beta = 0.0016
  const sb = 0.1
  const sd = 0.01

  const Utotal = [0.00001, 0.001, 0.0, 0.0, 0.00001]
  const msi_multiplier = 100.0
  const cin_multiplier = 0.01

  population = Dict{Array{Float64,1}, Lineage}()
  population[[0.0,0.0,0.0,1.0,0.0]] = Lineage(1.0,popN0,[0.0,0.0,0.0,1.0,0.0],1.0,0.0, 0)
    #initialize population - at the start everybody has fitness = 1.0, no mutations

  cancer_frequency = Float64
  cancer_frequency = cancer_state(population, 20)[1]

  msi_f, cin_f = assay_mutation_rate(population)
  popw, popN, varw = average_fitness(population)

  generations = Int64
  generations = 0

  cin_est_time = Int64
  cin_est_time = 0

  while cancer_frequency<.1
    popw, popN, varw = average_fitness(population)
 #   println("Generation: ", generations, " CIN frequency: ", cin_f, " Fitness: ", popw, " Population size: ", popN, " Cancer frequency: ", cancer_frequency)
    population = wright_fisher_reproduction(population, popN, round(Int,popN*(1+popw*beta)), popw)

    population = mutate_population(population, Utotal,  ben_dis, del_dis, sb, sd, msi_multiplier, cin_multiplier)

    cancer_frequency = cancer_state(population, 20)
    msi_f, cin_f = assay_mutation_rate(population)
    generations+=1
    if cin_f>=0.1
      if cin_est_time==0
        cin_est_time = generations
      end
    end
  end

  return generations, cin_est_time, cin_f
end
#end

function outfiles(time_to_mut, time_to_canc,cin_freq_canc)
  outfile = open(string(job_id,"time_to_mut.csv"), "a")
  write(outfile, join(time_to_mut, ","), "\n")
  close(outfile)
  outfile = open(string(job_id,"time_to_canc.csv"), "a")
  write(outfile, join(time_to_canc, ","), "\n")
  close(outfile)
  outfile = open(string(job_id,"cin_freq_canc.csv"), "a")
  write(outfile, join(cin_freq_canc, ","), "\n")
  close(outfile)
end


job_id = ARGS[1]

del_dis = 1
ben_dis = 4

time_to_mut= Int64[]
time_to_canc = Int64[]
cin_freq_canc = Float64[]
#
for run = 1:100
   tmp =simulate(ben_dis, del_dis)
   push!(time_to_mut, tmp[2])
   push!(time_to_canc, tmp[1])
   push!(cin_freq_canc, tmp[3])

end
outfiles(time_to_mut, time_to_canc,cin_freq_canc)
