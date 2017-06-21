arg = { ... }

assert(#arg == 4, "Expected 4 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = arg[1]
reverseOrbitTime = arg[1] / arg[2]

r0  = arg[3]

dwarfMass  = arg[4]

model1Bodies = 200001
totalBodies = model1Bodies
nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.32"

function makePotential()
  return nil
  -- return Potential.create{
  --     spherical = Spherical.spherical{ mass = 67479.9, scale = 0.6 },
  --     disk      = Disk.exponential{ mass = 224933, scaleLength = 4 },
  --     halo      = Halo.nfw{ vhalo = 120, scaleLength = 22.25 }
  -- }
end

encMass = plummerTimestepIntegral(r0, sqr(r0) + sqr(r0/arg[4]) , dwarfMass, 1e-5)

function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime * sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end
l = 218
b = 53.5
r = 28.6
-- Also required
function makeBodies(ctx, potential)
    local firstModel
    -- position = lbrToCartesian(ctx, Vector.create(l,b,r))
    position = Vector.create(0,0,0)
    local finalPosition, finalVelocity = position, Vector.create(0,0,0)

  -- local finalPosition, finalVelocity = reverseOrbit{
  --       potential = potential,
  --       position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
  --       velocity  = Vector.create(-156, 79, 107),
  --       tstop     = reverseOrbitTime,
  --       dt        = ctx.timestep / 10.0
  --   }
    firstModel = predefinedModels.plummer{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass        = dwarfMass,
        scaleRadius = r0,
        ignore      = true
    }

      for i,v in ipairs(firstModel)
      do
     	v.ignore = false
      end

return firstModel
end

function makeHistogram()
   return HistogramParams.create{
    phi = 128.79,
    theta = 54.39,
    psi = 90.70,
    lambdaStart = -10,
    lambdaEnd = 10,
    lambdaBins = 100,
    betaStart = -180,
    betaEnd = -180,
    betaBins = 1}
end


