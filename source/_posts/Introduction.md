---
title: Introduction
date: 2024-04-26 12:30:00
tags:
---
Welcome to my webpage! Here I document my progress and ideas as I develop my [Game](https://github.com/blackmagic919/CivGame/tree/main)!
I delve deeply into research and design choices regarding procedural content generation, and design architechture to optimize efficiency on [Unity's](https://unity.com/) platform.


## Basic Details

### Project Concept

![](Chunk_Borders.png)

[Procedural Generation](https://en.wikipedia.org/wiki/Procedural_generation) is a form of computer architechture that focuses on the dynamic generation of content. Contrary to traditional architechture, all broad-level states are not defined beforehand, but may be determined through a set of changing inputs. In terms of games, this often takes the form of terrain-generation, where near-infinite worlds may be simulated through the contained generation of only the viewer's immediate vicinity. This is a topic I will explore in-depth.

A common method to approach procedural generation is the tokenization of content which can be replaced selectively after its accessibility expires. Doing so reduces the load during generation of new content, as existing content may be cached while accessible. Regarding terrain-generation, tokenization is commonly referred to as chunks and divide continuous terrain into contained sections--a visualization of this concept can be seen above.

### Project Goals

While traditionally, 'games' require a storyline and conflict. I'm more content discovering and improving my implementations of procedural content generation. Doing so has allowed me to uncover incredible research surrounding data-management and GPU programming, and enabled me to advance my own skill-set through the chronic challenge entailed by its improvement.

Nevertheless, I do have ideas of eventually retrofitting this project in an actual play-able game. An interesting concept is [Aincrad from Sword Art Online](https://swordartonline.fandom.com/wiki/Aincrad), as generation is not limited in the vertical axis, a seemingly infinite amount of these 'levels' may be generated with increasing difficulty.

