arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)
print('seed=', argSeed)
evolveTime       = arg[1]
reverseOrbitTime = arg[1] / tonumber(arg[2])

r0  = arg[3]
light_r_ratio = arg[4]

dwarfMass  = arg[5]
light_mass_ratio = arg[6]

model1Bodies = 20000

totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.50"

-- 3.23e-4timestep
-- 2.25e-6soften

rscale_d = r0 / light_r_ratio
rscale_l = r0 
mass_d   = dwarfMass * (1.0 - light_mass_ratio)
mass_l   = dwarfMass * light_mass_ratio 


function makePotential()
   
return  nil
end

function get_timestep()
        --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3* ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    s1 = cube(rscale_l) / (mass_enc_d + mass_l)
    s2 = cube(rscale_d) / (mass_enc_l + mass_d)
    
    --return the smaller time step
    if(s1 < s2) then
        s = s1
    else
        s = s2
    end
    
    -- I did it this way so there was only one place to change the time step. 
    t = (1/100) * sqrt( pi_4_3 * s)
    print(t, t1, t2)
    return t
end





function makeContext()
   soften_length = (mass_l*rscale_l + mass_d*rscale_d)/(mass_d+mass_l)
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = get_timestep(),
      eps2       = calculateEps2(totalBodies, soften_length),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end


-- old verions of things, kept for posterity. 
-- -- encMass -- --
    -- mass_enc_l = mass_l* (rscale_l)^3* ( (sqr(rscale_l)+ sqr(rscale_l) ) )^(-3.0/2.0)
        -- calculates the mass of the light matter comp enclosed within its scale lenght. This is equally
        -- to, analytically, 1/(2sqrt(2)) of the total mass. 

    -- encMass = plummerTimestepIntegral(rscale_l,  rscale_d , mass_d, 1e-7)
        -- my fixed version. Mass of the dark comp enclosed within light matter scale r. 
        -- this is actually correct, but it calls a function when the value can be calculated directly with mass_enc_d

    -- encMass = plummerTimestepIntegral(r0*light_r_ratio, sqr(r0) + sqr(r0/light_r_ratio) , dwarfMass, 1e-7)
        -- this is incorrect we believe because of the input parameters. The function seems to try to calculate the mass of the dark 
        -- component that is enclosed within the light matter scale radius. This is not what is calculated with these inputs.

-- -- timestep and softening parameter: -- --
    --       timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
        -- this one is incorrect because it ovre calculates the denom by having the total dwarf mass and also the mass enclosed
        -- of the dark comp. 
    --       timestep   = (1/100.0) * sqrt((pi_4_3 * cube(rscale_l)) / (mass_enc_d + mass_enc_l)),
        -- this one is incorrect because it under calculates the denom. The mass enc of the light comp within the light scale r
        -- is a small fraction of the whole. 

    --       eps2       = calculateEps2(totalBodies, r0/light_r_ratio),
        -- this is incorrect because when we explore the parameter space and the light matter radius is small this becomes incorrect

-- Also required
function makeBodies(ctx, potential)
    local firstModel
    local finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)

    firstModel = predefinedModels.isotropic{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass1       = mass_l, 
        mass2       = mass_d,
        scaleRadius1 = rscale_l,
        scaleRadius2 = rscale_d,
        ignore      = true
    }

return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = 0.0,
     lambdaEnd = 300,
     lambdaBins = 200,
     betaStart = 0.0,
     betaEnd = -200,
     betaBins = 1
}
end

