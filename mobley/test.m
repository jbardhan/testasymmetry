function testset=test(i)

    ref_name=sprintf('OptWater_thermo_rand_%d',i);
    ref_data=load(ref_name);
    testset=ref_data.testset;
end
