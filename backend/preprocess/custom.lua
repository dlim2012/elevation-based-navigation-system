
local tables = {}

-- 'edges' table : table to preprocess OpenStreetMap "Way"
tables._edges = osm2pgsql.define_way_table('edges', {
    { column = 'type',     type = 'text', not_null = true},
    { column = 'oneway',   type = 'direction' },
    { column = 'nodes',    sql_type = 'int8[]'},
    { column = 'relation_ids',  sql_type = 'int8[]'},
    { column = 'geom',     projection = 4326, type = 'linestring', not_null = true },
    { column = 'tags',     type = 'jsonb' }, -- [todo]: remove (unused)
})

-- 'restriction' table: table to preprocess OpenStreetMap "Relation:restricion"
tables.restrictions = osm2pgsql.define_relation_table('restrictions', {
    { column = 'restriction', type = 'text' },
    { column = 'conditional', type = 'text' },
    { column = 'from_id', type = 'int8', not_null = true},
    { column = 'via_type', type='text', not_null = true},
    { column = 'via_id', type = 'int8', not_null = true },
    { column = 'to_id', type = 'int8', not_null = true },
    { column = 'tags',     type = 'jsonb' }, -- [todo]: remove (unused)
})

-- A data structure to save information for relation_ids when 'via' has a way type
-- OSM2PGSQL will _processNodes "Way" once more after populating w2r while processing "Relation"
local w2r = {}

-- Helper function to remove some of the tags we usually are not interested in.
-- Returns true if there are no tags left.
function clean_tags(tags)
    tags.odbl = nil
    tags.created_by = nil
    tags.source = nil
    tags['source:ref'] = nil

    return next(tags) == nil
end

-- highway types to save
-- local highway_types = {
--     'motorway',
--     'motorway_link',
--     'trunk',
--     'trunk_link',
--     'primary',
--     'primary_link',
--     'secondary',
--     'secondary_link',
--     'tertiary',
--     'tertiary_link',
--     'unclassified',
--     'residential',
--     'track',
--     'service',
--
--     'footway',
--     'path'
--     'steps',
--     'pedestrian',
--     'living_street',
--     'track',
--     'service',
--     'cycleway',
--     'road'
-- }

-- Prepare table "types" for quick checking of highway types
-- local types = {}
-- for _, k in ipairs(highway_types) do
--     types[k] = 1
-- end


-- Process "Ways"
function osm2pgsql.process_way(object)
--     -- We are only interested in highways
--     if not object.tags.highway then
--         return
--     end

    local row = {
        type = object.tags.highway or '',
        oneway = object.tags.oneway or 0,
        nodes = '{' .. table.concat(object.nodes, ',') .. '}',
        geom = object:as_linestring(),
        tags = object.tags
    }

    -- If there is any data from parent relations, add it in
    local d = w2r[object.id]
    if d then
        table.sort(d)
        row['relation_ids'] = '{' .. table.concat(d, ',') .. '}'
    end

    tables._edges:insert(row)
end

-- This function is called for every added, modified, or deleted relation.
-- Its only job is to return the ids of all member ways of the specified
-- relation we want to see in stage 2 again.
function osm2pgsql.select_relation_members(relation)
    -- Only interested in relations with type=route, route=road and a ref
    if relation.tags.type == 'restriction' then
        return { ways = osm2pgsql.way_member_ids(relation) }
    end
end

-- The process_relation() function should store all information about way
-- members that might be needed in stage 2.
function osm2pgsql.process_relation(object)
    -- Save only restriction
    if object.tags.type == 'restriction' then
        row = {
            id = object.id,
            restriction = object.tags.restriction,
            conditional = object.tags['restriction:conditional'],
            from_id = nil,
            to_id = nil,
            via_id = nil,
            via_type = nil,
            points = object.as_multilinestring(),
            tags = object.tags
        }

        -- Populate the members and prepare to add relation_ids in the ways table
        for _, member in ipairs(object.members) do
            if member.role == 'from' then
                row['from_id'] = member.ref
            elseif member.role == 'via' then
                row['via_id'] = member.ref
                if member.type ~= 'w' and member.type ~='n' then
                    return
                end
                if member.type == 'w' then
                    if not w2r[member.type] then
                        w2r[member.ref] = {}
                    end
                    w2r[member.ref][#(w2r[member.ref]) + 1] = object.id
                end
                row['via_type'] = member.type
            elseif member.role == 'to' then
                row['to_id'] = member.ref
            end

        end

        tables.restrictions:insert(row)
    end
end
